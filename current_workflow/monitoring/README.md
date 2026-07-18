# Monitoring Job

`refresh_live_report.sh` is the scheduled entrypoint for the live perturbation
report. It collects a fresh run snapshot, copies the canonical status and
history into this checkout, and invokes `report_generation_job.sh`.

`report_generation_job.sh` is the render-only entrypoint. It reads:

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


Run `refresh_live_report.sh` every 30 minutes from cron or a timer. Set
`PUBLISH_TO_GIT=1` to commit generated monitoring artifacts and push the
currently checked-out branch after each refresh. Source and configuration
files are excluded from automated commits. The cron example is in
`report_generation.cron`.
The history filename is legacy; the report job treats it as the canonical
monitor timeline.
