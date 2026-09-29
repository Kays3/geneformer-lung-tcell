"""Amendment 5 B1 fix: derive the 100-gene draw_order that every downstream script (run_phase8_null.sh,
build_ovx_index.py, ambient_loao_null.py, null_analysis.py) actually reads, from
frozen_100_estimable.json's `rows` (which has no draw_order key -- Stanley's finding). Asserts exactly
100 estimable rows and that the last row's position is 984, so a stale or wrong frozen file is refused
rather than silently producing the wrong N.
"""
import argparse
import json


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frozen", required=True, help="frozen_100_estimable.json")
    ap.add_argument("--expect-n", type=int, default=100)
    ap.add_argument("--expect-last-position", type=int, default=984)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    frozen = json.load(open(a.frozen))
    estimable = [r for r in frozen["rows"] if r["estimable"]]
    assert len(estimable) == a.expect_n, f"expected {a.expect_n} estimable genes, found {len(estimable)}"
    assert frozen["rows"][-1]["position"] == a.expect_last_position, (
        f"expected last row position {a.expect_last_position}, found {frozen['rows'][-1]['position']}")
    assert frozen["rows"][-1]["estimable"], "last row of the frozen prefix must itself be the 100th estimable gene"

    draw_order = [r["gene"] for r in estimable]
    out = {"about": __doc__.strip(), "source": a.frozen, "n": len(draw_order), "draw_order": draw_order}
    json.dump(out, open(a.out, "w"), indent=1)
    print(f"wrote {len(draw_order)} genes to {a.out}")


if __name__ == "__main__":
    main()
