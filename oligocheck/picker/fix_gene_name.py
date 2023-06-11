# %%
import re

import scanpy as sc

from oligocheck.boilerplate import jprint
from oligocheck.merfish.external_data import ExternalData

gtf = ExternalData(
    cache="data/mm39/gencode_vM32_transcripts.parquet",
    path="data/mm39/gencode.vM32.chr_patch_hapl_scaff.basic.annotation.gtf",
    fasta="data/mm39/combi.fa.gz",
)

adata = sc.read_h5ad("data/fpkm/embryo_raw.h5ad")
r = re.compile(r"Gm\d+")
converted, mapping, res = gtf.check_gene_names([x for x in adata.var_names if not r.match(x)])
sel = {}
# %%
dupes = {x[0] for x in res["dup"]}
for line in dupes:
    if line in sel:
        continue
    jprint(choices := {i: x for i, x in enumerate((x for x in res["out"] if x["query"] == line), 1)})
    inp = input()
    if not inp:
        break
    try:
        n = int(inp)
    except ValueError:
        sel[line] = inp
    else:
        if n == 0:
            sel[line] = line  # no idea
        else:
            sel[line] = choices[n]["symbol"]
print("Done!")

# %%
combi = sel | mapping
adata.var_names = [combi.get(x, x) for x in adata.var_names]

# %%
adata.write_h5ad("embryo_raw.h5ad")
# %%
# %%
