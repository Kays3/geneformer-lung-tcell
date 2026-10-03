# Graphical abstract

`graphical-abstract.png` (2400 × 1200) summarises whether the balanced-donor Geneformer design
transfers from lung (43 LUAD donors) to colon (E2, 19 Pelka 2021 donors), in three panels: per-donor
classifier balanced accuracy in both cohorts, the colon null-gene deletion/overexpression scatter
(H2b), and lung against colon deletion medians for the shared Panel B genes (H2c).

`build.py` reads only committed result files and checks every quoted number against the "Latest
studies" table in the top-level `README.md` before drawing. It recomputes the H2b Spearman rho from the
100 stored null-gene medians and the H2c sign agreement from the two `outcome_rows.json` files, and
stops if either differs from the registered result. The lung null-gene rho (−0.593) is quoted from the
LUAD report because its per-gene medians are not in the repository.

```sh
python3 build.py                            # writes graphical-abstract.svg and .html
NODE_PATH=$(npm root -g) node render.mjs    # writes graphical-abstract.png with Chromium (playwright)
```
