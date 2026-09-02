# PGCopilot benchmark: evaluation & analysis scripts

Reproducible evaluation pipeline for the manuscript
*PGCopilot: a proactive research agent for plant functional genomics*
(NPH-MET-2026-58029).

## Contents

| File | Purpose |
|---|---|
| `run_dify_evaluation.py` | Runs the multiple-choice benchmark through the Dify evaluation workflow and records per-question outputs for both modes (agent = LLM + MCP tools; LLM = no tools). |
| `eval_mcp_workflow_target.yml` | Exported Dify workflow DSL (target-plant app). Import into Dify, set the model-provider node to the LLM you want to evaluate. Contains **no credentials**. |
| `compute_summary_tables.py` | Recomputes all summary tables from the per-question output CSVs (Table 1 performance comparison with McNemar significance, species-stratified results with Wilson 95% CI, per-domain sample sizes). |
| `README.md` | This file. |

## Benchmark data

`../datasets/` (i.e. `workflow_evaluation/datasets/`) contains the
machine-readable benchmark (Expert MoBiPlant, Burda et al. 2025):

* `expert_mobi.json` – full 565-question benchmark
* `filtered_species.json` – 308-question target-plant subset
* `other_species.json` – 257-question non-target subset

Each entry: `question`, `options` (3 answer choices), `answer` (correct option
index), `plant_species`, `area` (research domain), `source` / `doi`.

## Running an evaluation

```bash
export DIFY_API_KEY=app-xxxxxxxx        # Dify app API key

python run_dify_evaluation.py \
    --benchmark ../datasets/filtered_species.json \
    --dataset target \
    --model-label DeepSeek-V3.2 \
    --output ../model_outputs/DeepSeek-V3.2_target.csv
```

One Dify app = one evaluated LLM. To evaluate a different model, import the
workflow YAML into Dify, switch the model-provider node, and run again.

Useful modes:

* `--mode all`      – run every question
* `--mode failed`   – resume: retry rows whose previous run did not succeed
* `--mode empty`    – resume: retry rows where the agent produced no answer

Results are written after every question, so interrupted runs can always be
resumed. Objective correctness (`agent_objective`, `llm_objective`) is computed
locally as `choice == correct_answer`.

## Reproducing the summary tables

```bash
python compute_summary_tables.py --inputs ../model_outputs/ --outdir ../summary_tables/
```

Produces:

* `Table1_Performance_Comparison.csv` – accuracy per model/dataset for both
  judgement paths (LLM-judge flag and objective choice-match), McNemar test
  (continuity-corrected) significance
* `Table_S_Species_Stratified.csv` – per-species accuracy with Wilson 95% CI
* `Table_S_Domain_SampleSizes.csv` – questions per research domain

## Data files

`model_outputs/` contains one CSV per model/dataset (6 models × 2 datasets).
Every row: question id, both answer modes, judge evaluation, objective
correctness, status and timing. A handful of rows have an empty agent answer
(see manuscript Methods: "persistently failing questions"); these are counted
as incorrect.
