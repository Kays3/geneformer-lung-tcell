# GSE263196 spatial validation

Tests whether the audit's pre-registered T-cell dysfunction signature is
enriched in T-cell-rich regions of SCLC tumor tissue, using the GSE263196
Visium cohort identified as an orthogonal spatial validation source by
`sclc_validation/audit/`.

## Result

Across 5 SCLC Visium samples (15,632 in-tissue spots total), a per-spot
T-cell marker score is positively correlated with a per-spot dysfunction
(exhaustion marker) score in 4 of 5 samples individually (Spearman
p < 1e-3), and in the pooled 5-sample meta-analysis: **ρ = 0.161, 95% CI
[0.146, 0.176]**. See `figures/tcell_dysfunction_correlation_forest.png`
and `METHODS.md` for the full design and scope decisions.

## Repository map

```text
spatial_validation.py   loads all 5 samples, scores spots, tests correlation
METHODS.md               design, scope decisions, limitations
results/                 per-sample and pooled correlation tables, spot-level scores
figures/                 forest plot of effect sizes with 95% CIs
```

Run:

```bash
/Users/kaisardauyey/workspace/research1/.venv/bin/python spatial_validation.py
```

Source Visium files remain outside Git (`../audit/source_metadata/GSE263196_RAW/`,
already downloaded and validated by the audit); compact result tables and the
figure are repository-trackable.
