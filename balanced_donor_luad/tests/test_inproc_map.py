"""inproc_map: with num_proc=1, Dataset.map/filter must run IN this process (the
Phase 6 slowdown was a forked one-worker pool per map), and give the same rows."""
import os
import sys

from datasets import Dataset

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import inproc_map  # noqa: E402


def _pids(ds):
    return set(ds.map(lambda x: {"pid": os.getpid()}, num_proc=1)["pid"])


def test_map_and_filter_run_in_process_after_install():
    inproc_map.install_inproc_map()
    ds = Dataset.from_dict({"a": list(range(20))})
    assert _pids(ds) == {os.getpid()}
    kept = ds.filter(lambda x: (x["a"] % 2 == 0) and os.getpid() == PARENT, num_proc=1)
    assert kept["a"] == list(range(0, 20, 2))


def test_unpatched_num_proc_1_forks_a_worker():
    # The failure mode: the original method with num_proc=1 runs in a different pid.
    orig = Dataset.map.__wrapped__ if hasattr(Dataset.map, "__wrapped__") else Dataset.map
    ds = Dataset.from_dict({"a": list(range(4))})
    pids = set(orig(ds, lambda x: {"pid": os.getpid()}, num_proc=1)["pid"])
    assert pids != {os.getpid()}


def test_install_is_idempotent():
    inproc_map.install_inproc_map(); inproc_map.install_inproc_map()
    assert Dataset.map.__wrapped__.__name__ == "map"


PARENT = os.getpid()
