# Immune axis test — is Normal < SCLC < LUAD a real ordering?

The poster and talk carry one unresolved claim: the model orders the three
states **Normal < SCLC < LUAD** on the checkpoint axis, which appears to
contradict the clinical picture of SCLC as an immune-cold, ICI-resistant tumour.
That claim rests on two shifts measured from a single source state.

This module tests it with the **overexpression arm**, which is census-complete
(every held-out source cell perturbed for every gene, N constant per state) and
therefore the only arm in which the three source states are directly comparable.

[**PLAN.md**](PLAN.md) — the full plan: what is being tested, what has been run,
the three analyses still needed, the replacement figure spec, and the results
that would overturn the conclusion.

## What runs today

```bash
python3 sclc_validation/immune_axis_test/axis_consistency.py
```

Laptop only — reads the committed stats tables from
[`../perturbation_workflow/targeted_panel/results/`](../perturbation_workflow/targeted_panel/results),
needs no GPU and no raw perturbation pickles. Writes six tables to `results/`.

Headline: **1 gene of 50 is consistent with a one-dimensional ordering**, and for
the seven-gene exhaustion program the "SCLC in the middle" prediction fails from
the Normal source for every gene. See PLAN.md §4 for the full result and §5 for
the working interpretation.

## Status

Test 1 run. Tests 2 (measured baseline expression), 3 (embedding geometry) and 4
(program-level overexpression, GPU) are specified but not executed. **No poster
or talk wording has been changed** — see PLAN.md §10.

```text
PLAN.md               the plan; read this first
axis_consistency.py   test 1, runs on the laptop
results/              its six output tables
```
