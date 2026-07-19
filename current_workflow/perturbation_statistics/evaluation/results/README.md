# Evaluation results workspace

Store future evaluation outputs here without modifying the completed primary
analysis tables.

Recommended layout:

```text
results/
├── analysis/
│   ├── tables/
│   └── figures/
└── biology/
    ├── tables/
    └── figures/
```

Every result set should record:

- source files and cohort version;
- generation timestamp and script or notebook;
- unit of analysis and donor handling;
- filters, thresholds, and multiple-testing method;
- model checkpoint and embedding layer;
- pass/fail outcome against the prespecified evaluation gate; and
- limitations that remain after the evaluation.
