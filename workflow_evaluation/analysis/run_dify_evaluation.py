#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_dify_evaluation.py
======================
Reproducible evaluation runner for the PGCopilot benchmark (Expert MoBiPlant subsets).

It sends multiple-choice plant-biology questions to a Dify workflow app
("eval mcp" evaluation workflow) and records, for every question:

  - with_mcp_answer / with_mcp_choice / with_mcp_evaluation / with_mcp_correct
    (agent mode: LLM + MCP tool access, i.e. PGCopilot)
  - no_mcp_answer  / no_mcp_choice  / no_mcp_evaluation  / no_mcp_correct
    (LLM mode: direct inference without tool augmentation)
  - objective correctness (choice == correct_answer), computed locally

The Dify app determines WHICH model is evaluated: import the exported workflow
DSL (`eval_mcp_workflow.yml`) into Dify and set the model provider node to the
desired LLM. Run this script once per model.

Workflow inputs (defined in the DSL):
  question : str   question text
  answer   : str   index (0/1/2) of the correct option
  options  : str   numbered option texts ("0. ...\n1. ...\n2. ...")

Benchmark files (Expert MoBiPlant, Burda et al. 2025):
  filtered_species.json : 308 target-plant questions
  other_species.json    : 257 non-target questions
  Each entry: {question, area, plant_species, options[3], answer(int), source, doi, ...}

Usage
-----
# 1. Full run (all questions of one dataset, one model per Dify app):
export DIFY_API_KEY=app-xxxxxxxx
python run_dify_evaluation.py \
    --benchmark ../datasets/filtered_species.json \
    --dataset target \
    --model-label DeepSeek-V3.2 \
    --output ../model_outputs/DeepSeek-V3.2_target.csv

# 2. Resume a crashed run, retrying failed rows only:
python run_dify_evaluation.py --output outputs/DeepSeek-V3.2_target.csv ...

# 3. Only rows whose agent answer is empty in an existing output file:
python run_dify_evaluation.py --mode empty --output outputs/Gemini_target.csv ...

Notes
-----
* Results are appended to --output after EVERY question (crash-safe resume).
* Never put your Dify API key in the code; use the environment variable.
"""
import argparse
import csv
import io
import json
import os
import re
import sys
import time

import requests

FIELDS = ['Model', 'Dataset', 'index', 'question', 'correct_answer',
          'with_mcp_answer', 'with_mcp_choice', 'with_mcp_evaluation', 'with_mcp_correct',
          'no_mcp_answer', 'no_mcp_choice', 'no_mcp_evaluation', 'no_mcp_correct',
          'status', 'elapsed_time', 'agent_objective', 'llm_objective']


def norm_text(s: str) -> str:
    """Normalise unicode spaces / artefacts so benchmark JSON matches question text."""
    for ch in ('\u2009', '\u00a0', '\u202f', '\ufffd'):
        s = s.replace(ch, ' ')
    return re.sub(r'\s+', ' ', s).strip()


def load_benchmark(path: str):
    with open(path, encoding='utf-8') as f:
        data = json.load(f)
    jmap = {}
    for e in data:
        jmap[norm_text(e['question'])] = e
    return jmap


def load_existing(path: str):
    if not os.path.exists(path):
        return [], set()
    with open(path, encoding='utf-8-sig', newline='') as f:
        rows = list(csv.DictReader(f))
    return rows, {(r['Dataset'], r['index']) for r in rows}


def blank_rec(**kw):
    rec = {c: '' for c in FIELDS}
    rec.update(kw)
    return rec


def objective(rec, side_prefix):
    ch = str(rec.get(side_prefix + '_choice') or '').strip()
    ca = str(rec.get('correct_answer') or '').strip()
    return 'True' if ch in ('0', '1', '2') and ch == ca else 'False'


def save_all(path, rows):
    with open(path, 'w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def run_one(session, base_url, api_key, payload, timeout):
    h = {'Authorization': 'Bearer ' + api_key, 'Content-Type': 'application/json'}
    resp = session.post(base_url.rstrip('/') + '/workflows/run',
                        json=payload, headers=h, timeout=timeout)
    d = resp.json().get('data', {})
    out = d.get('outputs', {}) or {}
    return d, out


def main():
    ap = argparse.ArgumentParser(description='PGCopilot benchmark evaluation via Dify workflow API')
    ap.add_argument('--benchmark', required=True, help='filtered_species.json / other_species.json')
    ap.add_argument('--dataset', required=True, choices=['target', 'nontarget'])
    ap.add_argument('--model-label', required=True, help='name of the evaluated LLM (recorded in output)')
    ap.add_argument('--output', required=True, help='output CSV (created / resumed)')
    ap.add_argument('--api-base', default=os.environ.get('DIFY_API_BASE', 'http://localhost/v1'))
    ap.add_argument('--api-key', default=os.environ.get('DIFY_API_KEY', ''),
                    help='Dify app API key (or set DIFY_API_KEY env var)')
    ap.add_argument('--mode', default='all', choices=['all', 'failed', 'empty'],
                    help='all=every question; failed=retry rows whose previous status != succeeded; '
                         'empty=only rows with empty with_mcp_answer (for "agent produced no answer" cases)')
    ap.add_argument('--limit', type=int, default=0, help='max questions this run (0 = no limit)')
    ap.add_argument('--timeout', type=int, default=300, help='per-request timeout (s)')
    ap.add_argument('--retries', type=int, default=1, help='retries per question on transport error')
    ap.add_argument('--delay', type=float, default=0.0, help='delay between questions (s)')
    args = ap.parse_args()

    if not args.api_key:
        sys.exit('ERROR: provide the Dify app API key via --api-key or DIFY_API_KEY')

    jmap = load_benchmark(args.benchmark)
    rows, done = load_existing(args.output)
    by_key = {(r['Dataset'], r['index']): r for r in rows}

    # decide which questions to (re)run (benchmark list order)
    todo = []
    with open(args.benchmark, encoding='utf-8') as f:
        blist = json.load(f)
    for i, e in enumerate(blist):
        key = (args.dataset, str(i))
        if key in done:
            prev = by_key.get(key)
            if args.mode == 'all':
                pass  # rerun everything
            elif args.mode == 'failed' and prev and prev['status'] == 'succeeded':
                continue
            elif args.mode == 'empty' and prev and str(prev.get('with_mcp_answer') or '').strip():
                continue
            elif args.mode == 'failed' or args.mode == 'empty':
                # remove previous failed record so the fresh one replaces it
                rows = [r for r in rows if (r['Dataset'], r['index']) != key]
        if args.limit and len(todo) >= args.limit:
            break
        todo.append((i, e))
    # 'all' mode must also respect limit
    if args.mode == 'all' and args.limit:
        todo = todo[:args.limit]

    print('Dataset=%s  Model=%s  mode=%s  questions to run: %d'
          % (args.dataset, args.model_label, args.mode, len(todo)), flush=True)

    session = requests.Session()
    n_ok = n_fail = 0
    for i, e in todo:
        q = e['question']
        ca = str(e['answer'])
        opts = '\n'.join('%d. %s' % (k, o) for k, o in enumerate(e['options']))
        rec = blank_rec(Model=args.model_label, Dataset=args.dataset, index=str(i),
                        question=q, correct_answer=ca)
        if norm_text(q) not in jmap and q not in jmap:
            rec['status'] = 'no_json'
        else:
            payload = {'inputs': {'question': q, 'answer': ca, 'options': opts},
                       'response_mode': 'blocking', 'user': 'benchmark-eval'}
            ok_run = False
            for attempt in range(args.retries + 1):
                try:
                    d, out = run_one(session, args.api_base, args.api_key, payload, args.timeout)
                    rec['status'] = d.get('status') or 'unknown'
                    rec['elapsed_time'] = d.get('elapsed_time') or 0
                    for k in ['with_mcp_answer', 'with_mcp_choice', 'with_mcp_evaluation',
                              'no_mcp_answer', 'no_mcp_choice', 'no_mcp_evaluation']:
                        rec[k] = out.get(k) or ''
                    for k in ['with_mcp_correct', 'no_mcp_correct']:
                        v = out.get(k)
                        rec[k] = '' if v is None else str(v)
                    ok_run = (rec['status'] == 'succeeded')
                    break
                except Exception as ex:
                    rec['status'] = 'api_error:' + str(ex)[:60]
                    if attempt < args.retries:
                        time.sleep(5)
        rec['agent_objective'] = objective(rec, 'with_mcp')
        rec['llm_objective'] = objective(rec, 'no_mcp')
        # replace any previous record for this key
        rows = [r for r in rows if (r['Dataset'], r['index']) != (args.dataset, str(i))]
        rows.append(rec)
        save_all(args.output, rows)
        done.add((args.dataset, str(i)))
        if rec['status'] == 'succeeded':
            n_ok += 1
        else:
            n_fail += 1
        print('[%s] idx %s  %s  agent_choice=%s agent_obj=%s' % (
            args.dataset, rec['index'], rec['status'], rec['with_mcp_choice'],
            rec['agent_objective']), flush=True)
        if args.delay:
            time.sleep(args.delay)

    print('\nDone. succeeded=%d  failed=%d  total_records=%d  -> %s'
          % (n_ok, n_fail, len(rows), args.output))


if __name__ == '__main__':
    main()
