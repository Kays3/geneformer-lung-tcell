# Relationship To KD

`KD` is the local data and generated-artifact workspace for this repository's
lung/T-cell Geneformer workflows. This repository contains version-controlled
methods, monitoring, validation, reporting, and migration tools; it reads
selected KD datasets, model checkpoints, statistics, and perturbation outputs.

The large files are intentionally outside Git in:

```text
/home/thinkstation2/workspace/KD
```

Keep the KD directory layout stable unless all workflow scripts and path
configuration are updated together.
