"""Budget watcher v3 for the resumed S100 LUAD run (Amendment 9): hard stop at
4.2 GPU-h TOTAL, counting the 0.202 h spent before the 04:14:20Z budget stop.

elapsed = prior spent (parsed from s100_luad_20260928.BUDGET_STOP_1, "at elapsed X h")
        + wall clock since the resume launch (s100_luad_20260928.started_utc, rewritten
          by the launcher at resume; the first launch time is kept in started_utc_1)
remaining = sum over delete/overexpress arms without a completion marker of
            37.88 s + 0.0394 s/cell (OLS on the 7 completed delete arms), rescaled by
            max(observed/projected on completed delete/overexpress arms, 1.0).
If elapsed + remaining > 4.2 h: SIGTERM the run's process group and write
s100_luad_20260928.BUDGET_STOP. Refuses to start unless <pgid> contains a
run_targeted_panel.py process. Usage: watch_s100_budget_v3.py <pgid>"""
import calendar, csv, json, os, re, signal, sys, time
from pathlib import Path
B = Path("/home/kaisar/workspace/geneformer-lung-tcell/sclc_validation/bf16_bench/runs")
MAN = B / "s100_luad_preflight_pergene_20260928/s100_luad_eligibility_manifest_pergene.csv"
LIMIT_S, A, PER = 4.2 * 3600, 37.88, 0.0394
pgid = int(sys.argv[1])
members = []
for d in Path("/proc").iterdir():
    if d.name.isdigit():
        try:
            if os.getpgid(int(d.name)) == pgid:
                members.append((d / "cmdline").read_bytes().replace(b"\0", b" ").decode(errors="replace"))
        except (ProcessLookupError, PermissionError, FileNotFoundError):
            pass
if not any("run_targeted_panel.py" in c for c in members):
    sys.exit(f"pgid {pgid} has no run_targeted_panel.py process: {members}")
prior_s = float(re.search(r"at elapsed ([\d.]+) h", (B / "s100_luad_20260928.BUDGET_STOP_1").read_text()).group(1)) * 3600
start = calendar.timegm(time.strptime((B / "s100_luad_20260928.started_utc").read_text().strip(), "%Y-%m-%dT%H:%M:%SZ"))
rows = [r for r in csv.DictReader(open(MAN)) if r["eligible"] == "True"]
cells = {r["gene"]: int(r["n_cells_after_cap"]) for r in rows}
arms = [(t, g) for t in ("delete", "overexpress") for g in cells]
proj = lambda g: A + PER * cells[g]
log = open(B / "s100_luad_20260928.budget_watch.log", "a")
log.write(f"{time.strftime('%FT%TZ', time.gmtime())} watcher v3 start: pgid={pgid} members={members} "
          f"prior_spent_h={prior_s/3600:.3f} resume_clock_from={time.strftime('%FT%TZ', time.gmtime(start))} limit_h=4.2\n"); log.flush()
while True:
    try:
        os.killpg(pgid, 0)
    except ProcessLookupError:
        log.write(f"{time.strftime('%FT%TZ', time.gmtime())} launcher gone; watcher v3 exits\n"); log.flush(); break
    obs = pj = remaining = 0.0; n_done = 0
    for t, g in arms:
        m = B / "s100_luad_20260928/targeted_panel/raw" / t / "luad" / f"targeted_luad_{g}.complete.json"
        if m.exists():
            obs += json.loads(m.read_text()).get("elapsed_seconds", 0.0); pj += proj(g); n_done += 1
        else:
            remaining += proj(g)
    ratio = obs / pj if pj else 1.0
    elapsed = prior_s + (time.time() - start)
    total = elapsed + remaining * max(ratio, 1.0)
    log.write(f"{time.strftime('%FT%TZ', time.gmtime())} v3 arms_done={n_done}/{len(arms)} elapsed_h={elapsed/3600:.3f} "
              f"rate_ratio={ratio:.3f} projected_total_h={total/3600:.3f}\n"); log.flush()
    if total > LIMIT_S:
        (B / "s100_luad_20260928.BUDGET_STOP").write_text(f"projected {total/3600:.3f} h > 4.2 at elapsed {elapsed/3600:.3f} h\n")
        os.killpg(pgid, signal.SIGTERM)
        log.write("BUDGET STOP (v3): sent SIGTERM to process group\n"); log.flush(); break
    time.sleep(60)
