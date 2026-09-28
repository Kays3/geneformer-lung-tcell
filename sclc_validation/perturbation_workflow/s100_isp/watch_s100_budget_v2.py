"""Budget watcher for the S100 LUAD run: hard stop at 3.4 GPU-h (Amendment 8).
Every 60 s: elapsed = wall clock since launch (conservative, includes loads);
remaining = projection for genes without a completion marker, rescaled by
observed/projected on finished arms. If elapsed + remaining > 3.4 h, kill the
launcher's process group and write a STOP file. Usage: watch_s100_budget.py <pgid>

v2 (2026-09-28): the clock starts at the launch time in
s100_luad_20260928.started_utc, not when the watcher starts, and the watcher
refuses to run unless <pgid> contains a run_targeted_panel.py process. v1 was
started with the wrong PGID (the ssh shell's group), so it could not have
stopped the run."""
import calendar, json, os, signal, sys, time
from pathlib import Path
import csv
B = Path("/home/kaisar/workspace/geneformer-lung-tcell/sclc_validation/bf16_bench/runs")
MAN = B / "s100_luad_preflight_pergene_20260928/s100_luad_eligibility_manifest_pergene.csv"
LIMIT_S, LOAD, PER = 3.4 * 3600, 6.7, (40.1 - 6.7) / 300
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
rows = [r for r in csv.DictReader(open(MAN)) if r["eligible"] == "True"]
cells = {r["gene"]: int(r["n_cells_after_cap"]) for r in rows}
panel = [r["gene"] for r in rows if r["role"] != "matched_control"]
arms = [("s100_luad_noop_20260928", "noop", g) for g in panel] + \
       [("s100_luad_20260928", t, g) for t in ("delete", "overexpress") for g in cells]
proj = lambda g: LOAD + PER * cells[g]
start = calendar.timegm(time.strptime((B / "s100_luad_20260928.started_utc").read_text().strip(), "%Y-%m-%dT%H:%M:%SZ"))
log = open(B / "s100_luad_20260928.budget_watch.log", "a")
log.write(f"{time.strftime('%FT%TZ', time.gmtime())} watcher v2 start: pgid={pgid} members={members} clock_from={time.strftime('%FT%TZ', time.gmtime(start))}\n"); log.flush()
while True:
    try:
        os.killpg(pgid, 0)
    except ProcessLookupError:
        log.write(f"{time.strftime('%FT%TZ', time.gmtime())} launcher gone; watcher exits\n"); log.flush(); break
    done_obs = done_proj = 0.0; remaining = 0.0; n_done = 0
    for tag, t, g in arms:
        m = B / tag / "targeted_panel/raw" / t / "luad" / f"targeted_luad_{g}.complete.json"
        if m.exists():
            e = json.loads(m.read_text()).get("elapsed_seconds", 0.0)
            done_obs += e; done_proj += proj(g); n_done += 1
        else:
            remaining += proj(g)
    ratio = (done_obs / done_proj) if done_proj else 1.0
    elapsed = time.time() - start
    total = elapsed + remaining * max(ratio, 1.0)
    log.write(f"{time.strftime('%FT%TZ', time.gmtime())} arms_done={n_done}/{len(arms)} elapsed_h={elapsed/3600:.3f} "
              f"rate_ratio={ratio:.3f} projected_total_h={total/3600:.3f}\n"); log.flush()
    if total > LIMIT_S:
        (B / "s100_luad_20260928.BUDGET_STOP").write_text(f"projected {total/3600:.3f} h > 3.4 at elapsed {elapsed/3600:.3f} h\n")
        os.killpg(pgid, signal.SIGTERM)
        log.write("BUDGET STOP: sent SIGTERM to process group\n"); log.flush(); break
    time.sleep(60)
