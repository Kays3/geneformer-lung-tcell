# S100 LUAD analysis provenance (2026-09-28)

- Stanley's post-run PASS: ~17:30 JST. Michael's GO: 17:27:52 JST.
- Analysis start: 17:30:23 JST (thinkstation1, creation of runs/s100_luad_analysis_20260928/). Outputs written 17:30:45-17:30:53 JST.
- The four stats CSVs (delete/overexpress x luad_to_{normal,sclc}; mtimes 17:22:37-17:22:41 JST, written by the runner) were first read by this analysis: atime 17:30:44 JST. Their atime was <= mtime at Stanley's check.
- The two no-op stats CSVs were not read (atime 13:07:54, still before mtime); no-op scores come from the raw pickles.
- Code: analyse_s100_luad_20260928.py plus the analysis layer at 295aa2c (on top of chain head d7982063). The one code change after the chain head is the find_raw_pickle no-op glob fix and noop_raw_root; see that commit.
- Inputs: runner 974535b7, run panel 8a0b3668, control sets 607ac6c9; raw roots s100_luad_20260928 and s100_luad_noop_20260928.
