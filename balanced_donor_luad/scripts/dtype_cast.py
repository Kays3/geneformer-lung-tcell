"""Shared bf16/fp32 dtype-cast wrapper for the ISP runner scripts.

InSilicoPerturber (geneformer/in_silico_perturber.py) does
`from . import perturber_utils as pu` and calls `pu.load_model(...)`
internally -- our runners never call load_model directly, so there is no
call site to add a --dtype argument to without touching upstream geneformer
source (out of scope; see geneformer-bf16-export.patch for the one upstream
change we do make, which is unrelated: the 4-site .float() numpy-export
fix). Instead we monkeypatch the `load_model` attribute on the
`perturber_utils` module object itself -- `pu` is that same module object,
so the patch is visible through the alias with no signature change anywhere
upstream.

install_dtype_cast() must run before the first perturb_data() call in a
process (it patches at call time, not import time, so it does not need to
run before `from geneformer import InSilicoPerturber`).
"""
from __future__ import annotations

DTYPES = ("fp32", "bf16")


def install_dtype_cast(dtype: str) -> None:
    if dtype not in DTYPES:
        raise ValueError(f"dtype must be one of {DTYPES}, got {dtype!r}")

    import torch
    from geneformer import perturber_utils as pu

    if getattr(pu.load_model, "_bf16_bench_wrapped", False):
        return  # already installed in this process (e.g. re-entrant call)

    original_load_model = pu.load_model

    def load_model_cast(*args, **kwargs):
        model = original_load_model(*args, **kwargs)
        if dtype == "bf16":
            model = model.to(torch.bfloat16)
        return model

    load_model_cast._bf16_bench_wrapped = True
    load_model_cast._bf16_bench_dtype = dtype
    pu.load_model = load_model_cast
