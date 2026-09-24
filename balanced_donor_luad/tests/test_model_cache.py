"""model_cache's store must pickle as EMPTY, so that datasets' dill fingerprinting of
every .map() does not serialise the cached models (the Phase 6 slowdown, 2026-09-24).
A plain dict fails this test by construction."""
import os
import pickle
import sys

import dill

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import model_cache as mc  # noqa: E402

BIG = bytes(50_000_000)


def test_store_pickles_empty_even_when_full():
    s = mc._OpaqueStore(); s["k"] = BIG
    assert len(pickle.dumps(s)) < 1000 and len(dill.dumps(s)) < 1000


def test_store_behaves_as_dict_in_process():
    s = mc._OpaqueStore(); s["k"] = 1
    assert s["k"] == 1 and "k" in s and pickle.loads(pickle.dumps(s)) == {}


def test_plain_dict_would_fail():
    d = {"k": BIG}
    assert len(pickle.dumps(d)) > 1000   # the failure mode the store exists to prevent


def test_module_cache_is_the_opaque_store():
    assert isinstance(mc._CACHE, mc._OpaqueStore)
