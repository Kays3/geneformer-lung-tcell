# Graphical abstract

`graphical-abstract.png` (2400 × 1200) asks whether the evaluation criteria of the Geneformer protocol
(Zhang, Venkatesh & Theodoris, *Nature Protocols* 2026) hold in T cells under donor-level tests, with one
panel per criterion: held-out-patient classification (per-donor balanced accuracy in 43 lung and 19 colon
donors), the comparison with random genes (the colon null-gene deletion/overexpression scatter, H2b), and
candidate genes (lung against colon deletion medians for the shared Panel B genes, H2c). The footer points
to E3, the designed study of the same criteria in T cells of the protocol's own heart atlas.

`build.py` reads only committed result files and checks every quoted number against the "Latest
studies" table in the top-level `README.md` before drawing. It recomputes the H2b Spearman rho from the
100 stored null-gene medians and the H2c sign agreement from the two `outcome_rows.json` files, and
stops if either differs from the registered result. The lung null-gene rho (−0.593) is quoted from the
LUAD report because its per-gene medians are not in the repository.

```sh
python3 build.py                            # writes graphical-abstract.svg and .html
NODE_PATH=$(npm root -g) node render.mjs    # writes graphical-abstract.png with Chromium (playwright)
```
