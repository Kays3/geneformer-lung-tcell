# Monitoring Job

`report_generation_job.sh` is the entrypoint for the live perturbation report.
It reads:

- `latest_status.json`
- `hourly_history.csv`

and rewrites:

- `GPU_PROGRESS_REPORT.md`
- `progress_animation.gif`
- `progress_animation.svg`
- `cell_interaction_diagram.svg`
- `snapshot_gallery/*.svg`

The report also inspects the six expected directional CSV files in
`PERTURBATION_STATS_DIR` and reports each comparison as waiting, queued/running,
complete, empty, or unreadable. `report_generation_job.sh` supplies the live
run's statistics directory by default; set the environment variable to use a
different run.


Run it every 15 minutes from cron, a timer, or the existing workflow wrapper.
The cron example is in `report_generation.cron`.
The history filename is legacy; the report job treats it as the canonical
monitor timeline.
