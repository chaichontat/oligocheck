# %%
from pathlib import Path

import numpy as np
import polars as pl

from oligocheck.external.external_data import ExternalData
from oligocheck.picker.codebook import CodebookPicker

gtf_all = ExternalData(
    cache="data/mm39/gencode_vM32_transcripts_all.parquet",
    path="data/mm39/Mus_musculus.GRCm39.109.gtf",
    fasta="data/mm39/combi.fa.gz",
)

gtf = ExternalData(
    cache="data/mm39/gencode_vM32_transcripts.parquet",
    path="data/mm39/gencode.vM32.chr_patch_hapl_scaff.basic.annotation.gtf",
    fasta="data/mm39/combi.fa.gz",
)
# %%

pl.Config.set_fmt_str_lengths(100)
pl.Config.set_tbl_rows(100)

dfs, sams, filtereds, offtargets, overlapped = {}, {}, {}, {}, {}
for file in Path("output").glob("*_all.parquet"):
    gene = file.stem.split("_")[0]
    dfs[gene] = pl.read_parquet(f"output/{gene}_final.parquet")
    sams[gene] = pl.read_parquet(f"output/{gene}_all.parquet")
    filtereds[gene] = pl.read_parquet(f"output/{gene}_filtered.parquet")
    offtargets[gene] = pl.read_parquet(f"output/{gene}_offtargets.parquet")
    try:
        overlapped[gene] = pl.read_parquet(f"output/{gene}_final_overlapped.parquet")
    except FileNotFoundError:
        pass

mhd4 = CodebookPicker(np.loadtxt("data/readouts/mhd4_22bit.csv", delimiter=",", dtype=bool))
fpkm = pl.read_parquet("data/fpkm/P0_combi.parquet")

# %%
transcripts = (
    pl.DataFrame({"gene": dfs.keys()})
    .with_columns(pl.col("gene").apply(gtf.get_transcripts).alias("transcript"))
    .explode("transcript")
    .join(fpkm, left_on="transcript", right_on="transcript_id(s)", how="left")
    .join(
        gtf_all[["transcript_id", "transcript_name"]],
        left_on="transcript",
        right_on="transcript_id",
    )
    .fill_null(0)
    .groupby("gene")
    .agg(pl.col("FPKM").sum().alias("FPKM"))
)
# %%
# %%
import scanpy as sc

singlecell = sc.read_h5ad("data/fpkm/BrainAgingSpatialAtlas_snRNAseq.h5ad")
percentile = np.percentile(singlecell.X, 99.9, axis=0)
# %%
counts = dict(
    zip(singlecell.var_names, np.clip(singlecell.var["mean"] + percentile / singlecell.var["std"], 0, np.inf))
)
# %%
# singlecell = dict(pl.read_csv("data/fpkm/combi99percentile.csv").iter_rows())
# %%

counts = {x: counts.get(gtf.gene_to_eid(x), 0) for x in dfs.keys()}


# %%
def get_canonical(gene: str):
    return gtf_all.filter(pl.col("transcript_name") == f"{gene}-201")[0, "transcript_id"]


# counts = (
#     pl.DataFrame({"transcript": [get_canonical(x) for x in dfs.keys()]})
#     .join(fpkm, left_on="transcript", right_on="transcript_id(s)", how="left")
#     .join(
#         gtf_all[["transcript_id", "transcript_name"]],
#         left_on="transcript",
#         right_on="transcript_id",
#     )
# )
# %%

best, values = mhd4.find_optimalish(np.array(list(counts.values()))[:, np.newaxis], iterations=10000)
# %%
codebook = mhd4.gen_codebook(best).astype(np.uint8)
# pl.DataFrame(codebook).with_columns(gene=list(counts.keys()))
df = (
    pl.DataFrame(codebook)
    .rename({f"column_{i}": f"bit{i+1}" for i in range(codebook.shape[1])})
    .with_columns(
        gene=pl.Series(
            list(counts.keys()) + [f"Blank-{i+1}" for i in range(codebook.shape[0] - len(counts))]
        ),
    )
)
df = df.select("gene", *df.columns[:-1])
df.write_csv("panels/celltype_codebook.csv")


# %%

# %%


# %%
