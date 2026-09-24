"""Load each fine-tuned model once per process for the ~15,000 per-donor ISP calls.

InSilicoPerturber and EmbExtractor both obtain their model through
`perturber_utils.load_model(model_type, num_classes, model_directory, mode, ...)`.
This wrapper memoises that call on its arguments, in the same monkeypatch
style as bf16_bench/dtype_cast.py (install dtype_cast FIRST, then this, so
the cached object is the already-cast bf16 model). It does not touch any
Geneformer arithmetic: the cached object is exactly what the original
function returned on its first call.

Safety: ISP and embedding extraction only run inference. The wrapper still
records the model's parameter checksum at first load and `verify()` recomputes
it, so any in-place mutation during a run is detected rather than assumed away.
"""
from __future__ import annotations

import hashlib



class _OpaqueStore(dict):
    """A dict that pickles as EMPTY. Hugging Face `datasets` fingerprints every .map() by
    dill-pickling the mapped function and what it reaches; that walk reached this cache and
    serialised every cached 316M model on every map (~40-60 s per ISP call, CPU-bound).
    Pickling the store as empty removes that cost. It affects only the fingerprint, never a
    value computed from a model. Each donor's ISP input only ever meets that donor's one fold
    model, so no datasets cache entry can be shared across models."""
    def __reduce__(self):
        return (_OpaqueStore, ())


_CACHE: dict = _OpaqueStore()
_CHECKSUM: dict = {}


def _param_checksum(model) -> str:
    h = hashlib.sha256()
    for name, p in sorted(model.state_dict().items()):
        h.update(name.encode())
        h.update(p.detach().float().cpu().numpy().tobytes())
    return h.hexdigest()


def install_model_cache(checksum: bool = True) -> None:
    from geneformer import perturber_utils as pu

    if getattr(pu.load_model, "_model_cache_wrapped", False):
        return
    original = pu.load_model

    def cached_load_model(*args, **kwargs):
        key = (args, tuple(sorted(kwargs.items())))
        if key not in _CACHE:
            model = original(*args, **kwargs)
            _CACHE[key] = model
            if checksum:
                _CHECKSUM[key] = _param_checksum(model)
        return _CACHE[key]

    cached_load_model._model_cache_wrapped = True
    pu.load_model = cached_load_model


def verify() -> dict:
    """Recompute every cached model's checksum; raise if any changed."""
    changed = [k for k, m in _CACHE.items() if k in _CHECKSUM and _param_checksum(m) != _CHECKSUM[k]]
    if changed:
        raise RuntimeError(f"cached model(s) mutated during the run: {changed}")
    return {"n_models": len(_CACHE), "verified": True}


def cache_size() -> int:
    return len(_CACHE)
