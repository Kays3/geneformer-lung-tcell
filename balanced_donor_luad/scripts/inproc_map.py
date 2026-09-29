"""Run Hugging Face `datasets` map/filter IN PROCESS when num_proc == 1.

This `datasets` version treats num_proc=1 as "a pool of one worker". Every
Dataset.map/filter inside Geneformer's InSilicoPerturber (called with nproc=1)
then forks the whole ISP process (~20 GB resident: five cached 316M models) and
pipes results back through a multiprocess manager: ~6 s per map, several maps
per call. num_proc=None is the documented in-process path. The mapped function,
its inputs and its outputs are unchanged; only where it executes changes.
Install once per process, after dtype_cast and model_cache.
"""
from __future__ import annotations


def _wrap(method):
    def wrapper(self, *args, **kwargs):
        if kwargs.get("num_proc") == 1:
            kwargs["num_proc"] = None
        return method(self, *args, **kwargs)
    wrapper._inproc_wrapped = True
    wrapper.__wrapped__ = method
    return wrapper


def install_inproc_map() -> None:
    from datasets import Dataset
    for name in ("map", "filter"):
        m = getattr(Dataset, name)
        if not getattr(m, "_inproc_wrapped", False):
            setattr(Dataset, name, _wrap(m))
