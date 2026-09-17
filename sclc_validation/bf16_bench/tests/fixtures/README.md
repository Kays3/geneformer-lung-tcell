Verbatim `geneformer/emb_extractor.py` and `geneformer/perturber_utils.py`
from the pinned checkout (`f45a6c7de57ff07f946f146c254da02a90e2cdf5`,
`/home/kaisar/workspace/geneformer-uv-starter/Geneformer` on
thinkstation1), fetched read-only for `test_patch_idempotence.py`'s
fixture. Not modified. Confirmed identical to what both runner scripts
actually import at runtime (`/srv/lab/geneformer`, currently 4 commits
ahead at `04c2b2e` -- diffed and these two files are byte-for-byte the
same at both commits; see the bf16-replication task's Step 0 report for
the full checkout-pin discussion).
