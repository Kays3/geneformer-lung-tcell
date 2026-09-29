# Phase 6 run notes (ISP, both hosts)

Written 2026-09-25T00:22:37Z. All times are file times or `date -u` reads.

## Segments
Each resume appends one line to `phase6/started_utc_<host>.txt` and `phase6/code_commit_<host>.txt`; nothing is overwritten.

| Segment | Host | Start | Stop | Code | Calls done at stop |
|---|---|---|---|---|---|
| 1 | ts1 | 23:25:53Z | 23:49:15Z | fb2a6ef | 24 |
| 1 | ts2 | 23:25:55Z | 23:49:24Z | fb2a6ef | 24 |
| 2 | ts1 | 23:50:14Z | 00:13:25Z | 27155c3 | 120 |
| 2 | ts2 | 23:50:17Z | 00:13:26Z | 27155c3 | 112 |
| 3 | both | ~00:14Z (see started_utc) | — | 2cfe3eb | — |

Segments 1 and 2 were stopped deliberately, to apply two CPU-side performance fixes:
- 27155c3: the model cache pickles as empty for `datasets` fingerprinting;
- 2cfe3eb: `datasets` map/filter runs in process when num_proc=1.

Neither fix touches arithmetic. Evidence:
- Test calls: 8/8 and then 20/20 output pickles byte-identical to the pre-fix code.
- **Real stored outputs:** 12 stored pickles per host were recomputed under 2cfe3eb in a separate folder, and all 24 match by sha256. By writing segment: seg1 4+4, seg2 7+7, seg3 1+1.

## Reading the run log
- A resumed run logs a second `gene_done` line for genes it skipped entirely. Count completed work from the `<op>/<gene>/<donor>.complete.json` markers, never from log lines.
- Each marker is written after its pickle (run_isp.py: perturb_data, then json.dump of the marker), so a partial pickle never has a marker.
- The marker write itself is not atomic. The analysis rejects any marker that does not parse as JSON, and treats that call as not run.

## Failure policy
One failing call halts that host (`set -e`, no per-call try/except). This is deliberate: a halt is visible and resumable, whereas skip-and-record would leave silent holes. A (gene, donor) pair that fails deterministically goes to god as a decision.
