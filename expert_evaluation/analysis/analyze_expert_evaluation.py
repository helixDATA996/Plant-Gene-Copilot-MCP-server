#!/usr/bin/env python
"""
Reproducible analysis of the expert evaluation.

Default input:
    ../data/expert_scores_anonymized.xlsx

Outputs:
    expert_statistics.xlsx
    inter_rater_reliability.csv
    pairwise_significance_tests.csv

The script calculates:
    - Overall and stratified model performance
    - ICC(2,1) and ICC(2,k)
    - Paired Wilcoxon signed-rank tests
    - Paired t-tests
    - Cohen's dz

Required packages:
    pandas, numpy, scipy, openpyxl
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import ttest_rel, wilcoxon


MODEL_COLUMNS = [
    "Full PGC Agent",
    "PlantGPT",
    "Ablated General LLM",
    "PGC Agent without System Prompt",
]


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    default_input = script_dir.parent / "data" / "expert_scores_anonymized.xlsx"
    default_output = script_dir

    parser = argparse.ArgumentParser(
        description="Calculate reproducible expert-evaluation statistics."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=default_input,
        help="Input expert-score Excel file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=default_output,
        help="Directory for generated result files.",
    )
    return parser.parse_args()


def load_scores(path: Path) -> pd.DataFrame:
    """Load the long-format Reviewed Scores sheet and validate its schema."""
    workbook = pd.ExcelFile(path)
    sheet = "Reviewed Scores" if "Reviewed Scores" in workbook.sheet_names else workbook.sheet_names[0]
    df = pd.read_excel(path, sheet_name=sheet)

    required = {"Evaluator", "Question", "Metric", *MODEL_COLUMNS}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(
            f"Missing required columns: {sorted(missing)}. "
            f"Available columns: {list(df.columns)}"
        )

    df = df.dropna(subset=["Evaluator", "Question", "Metric"]).copy()
    for model in MODEL_COLUMNS:
        df[model] = pd.to_numeric(df[model], errors="raise")
        if ((df[model] < 0) | (df[model] > 10)).any():
            raise ValueError(f"Scores for {model!r} must be between 0 and 10.")

    df["Evaluator"] = df["Evaluator"].astype(str)
    df["Question"] = df["Question"].astype(str).str.lower()
    df["Metric"] = df["Metric"].astype(str).str.upper()

    key = ["Evaluator", "Question", "Metric"]
    duplicates = df.duplicated(key, keep=False)
    if duplicates.any():
        examples = df.loc[duplicates, key].head().to_dict("records")
        raise ValueError(
            "Duplicate evaluator-question-metric records detected. "
            f"Examples: {examples}"
        )
    return df


def population_sd(values: Iterable[float]) -> float:
    values = np.asarray(list(values), dtype=float)
    return float(np.std(values, ddof=0))


def summary_for(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for model in MODEL_COLUMNS:
        values = df[model].to_numpy(dtype=float)
        rows.append(
            {
                "Model": model,
                "N": len(values),
                "Mean score": round(float(np.mean(values)), 3),
                "SD": round(population_sd(values), 3),
                "Minimum": float(np.min(values)),
                "Maximum": float(np.max(values)),
            }
        )
    result = pd.DataFrame(rows).sort_values(
        "Mean score", ascending=False, ignore_index=True
    )
    result.insert(0, "Rank", np.arange(1, len(result) + 1))
    return result


def icc_two_way_random(matrix: np.ndarray) -> tuple[float, float]:
    """
    Calculate ICC(2,1) and ICC(2,k) using Shrout-Fleiss notation.

    Rows are targets and columns are raters. This is the two-way
    random-effects, absolute-agreement model.
    """
    x = np.asarray(matrix, dtype=float)
    n, k = x.shape
    if n < 2 or k < 2:
        return math.nan, math.nan

    grand_mean = x.mean()
    row_means = x.mean(axis=1)
    col_means = x.mean(axis=0)

    ss_rows = k * np.sum((row_means - grand_mean) ** 2)
    ss_cols = n * np.sum((col_means - grand_mean) ** 2)
    ss_total = np.sum((x - grand_mean) ** 2)
    ss_error = ss_total - ss_rows - ss_cols

    ms_rows = ss_rows / (n - 1)
    ms_cols = ss_cols / (k - 1)
    ms_error = ss_error / ((n - 1) * (k - 1))

    icc21 = (ms_rows - ms_error) / (
        ms_rows
        + (k - 1) * ms_error
        + k * (ms_cols - ms_error) / n
    )
    icc2k = (ms_rows - ms_error) / (ms_rows + (ms_cols - ms_error) / n)
    return float(icc21), float(icc2k)


def reliability_table(df: pd.DataFrame) -> pd.DataFrame:
    """Calculate overall, by-metric, and by-question inter-rater reliability."""
    long = df.melt(
        id_vars=["Evaluator", "Question", "Metric"],
        value_vars=MODEL_COLUMNS,
        var_name="Model",
        value_name="Score",
    )
    long["Target"] = (
        long["Question"] + "_" + long["Metric"] + "_" + long["Model"]
    )

    records = []

    def add_record(scope: str, question: str, metric: str, subset: pd.DataFrame) -> None:
        matrix = subset.pivot_table(
            index="Target", columns="Evaluator", values="Score", aggfunc="mean"
        ).dropna(axis=0)
        icc21, icc2k = icc_two_way_random(matrix.to_numpy())
        records.append(
            {
                "Scope": scope,
                "Question": question,
                "Metric": metric,
                "Targets": matrix.shape[0],
                "Raters": matrix.shape[1],
                "ICC(2,1)": round(icc21, 3),
                "ICC(2,k)": round(icc2k, 3),
            }
        )

    add_record("Overall", "All", "All", long)
    for metric in sorted(long["Metric"].unique(), key=lambda x: int(x[1:])):
        add_record("By metric", "All", metric, long[long["Metric"] == metric])
    for question in sorted(long["Question"].unique(), key=lambda x: int(x[1:])):
        add_record(
            "By question",
            question.upper(),
            "All",
            long[long["Question"] == question],
        )
    return pd.DataFrame(records)


def significance_table(df: pd.DataFrame) -> pd.DataFrame:
    """Compare models on identical evaluator-question-metric records."""
    records = []
    reference = MODEL_COLUMNS[0]

    for comparison in MODEL_COLUMNS[1:]:
        x = df[reference].to_numpy(dtype=float)
        y = df[comparison].to_numpy(dtype=float)
        differences = x - y

        wilcoxon_result = wilcoxon(
            x, y, zero_method="wilcox", alternative="two-sided"
        )
        t_result = ttest_rel(x, y)
        difference_sd = np.std(differences, ddof=1)
        cohens_dz = np.mean(differences) / difference_sd

        records.append(
            {
                "Model A": reference,
                "Model B": comparison,
                "N paired ratings": len(differences),
                "Mean A": round(float(np.mean(x)), 3),
                "Mean B": round(float(np.mean(y)), 3),
                "Mean difference A-B": round(float(np.mean(differences)), 3),
                "Wilcoxon statistic": round(float(wilcoxon_result.statistic), 3),
                "Wilcoxon p": float(wilcoxon_result.pvalue),
                "Paired t statistic": round(float(t_result.statistic), 3),
                "Paired t-test p": float(t_result.pvalue),
                "Cohen's dz": round(float(cohens_dz), 3),
                "A higher than B": int(np.sum(differences > 0)),
                "Tie": int(np.sum(differences == 0)),
                "A lower than B": int(np.sum(differences < 0)),
            }
        )
    return pd.DataFrame(records)


def write_results(
    output_dir: Path,
    input_records: int,
    overall: pd.DataFrame,
    by_question: pd.DataFrame,
    by_metric: pd.DataFrame,
    reliability: pd.DataFrame,
    significance: pd.DataFrame,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    workbook_path = output_dir / "expert_statistics.xlsx"

    with pd.ExcelWriter(workbook_path, engine="openpyxl") as writer:
        overall.to_excel(writer, sheet_name="Overall Performance", index=False)
        by_question.to_excel(writer, sheet_name="Performance by Question", index=False)
        by_metric.to_excel(writer, sheet_name="Performance by Metric", index=False)
        reliability.to_excel(writer, sheet_name="Inter-rater Reliability", index=False)
        significance.to_excel(writer, sheet_name="Significance Tests", index=False)

        notes = pd.DataFrame(
            [
                ["Input file", "Expert score workbook, Reviewed Scores sheet"],
                ["ICC model", "Two-way random-effects, absolute agreement: ICC(2,1) and ICC(2,k)"],
                ["Significance test", "Paired Wilcoxon signed-rank test; paired t-test as a reference"],
                ["Effect size", "Cohen's dz = mean paired difference / SD of paired differences"],
                ["Model order", "; ".join(MODEL_COLUMNS)],
            ],
            columns=["Item", "Description"],
        )
        notes.to_excel(writer, sheet_name="Statistics Notes", index=False)

    reliability.to_csv(
        output_dir / "inter_rater_reliability.csv",
        index=False,
        encoding="utf-8-sig",
    )
    significance.to_csv(
        output_dir / "pairwise_significance_tests.csv",
        index=False,
        encoding="utf-8-sig",
    )

    print(f"Input records: {input_records}")
    print(f"Created: {workbook_path}")
    print(f"Created: {output_dir / 'inter_rater_reliability.csv'}")
    print(f"Created: {output_dir / 'pairwise_significance_tests.csv'}")


def main() -> None:
    args = parse_args()
    scores = load_scores(args.input)

    overall = summary_for(scores)
    by_question = (
        scores.groupby("Question", sort=False)
        .apply(summary_for, include_groups=False)
        .reset_index(level=1, drop=True)
        .reset_index()
    )
    by_metric = (
        scores.groupby("Metric", sort=False)
        .apply(summary_for, include_groups=False)
        .reset_index(level=1, drop=True)
        .reset_index()
    )

    # The grouped summaries are useful internally, but flatten their labels.
    by_question = by_question.rename(columns={"Question": "Question"})
    by_metric = by_metric.rename(columns={"Metric": "Metric"})

    reliability = reliability_table(scores)
    significance = significance_table(scores)
    write_results(
        args.output_dir,
        len(scores),
        overall,
        by_question,
        by_metric,
        reliability,
        significance,
    )


if __name__ == "__main__":
    main()
