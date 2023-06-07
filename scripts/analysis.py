# %%
import logging
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import polars as pl

from oligocheck.merfish.external_data import ExternalData
from oligocheck.merfish.nnupack import gen_model, nonspecific_test, secondary_structure
from oligocheck.seqcalc import tm_fish, tm_hybrid, tm_match
from oligocheck.sequtils import reverse_complement, stripplot

gtf_all = ExternalData(
    cache="data/mm39/gencode_vM32_transcripts_all.parquet",
    path="data/mm39/Mus_musculus.GRCm39.109.gtf",
    fasta="data/mm39/combi.fa.gz",
)
pl.Config.set_fmt_str_lengths(100)
pl.Config.set_tbl_rows(100)

dfs, sams, filtereds, offtargets, overlapped = {}, {}, {}, {}, {}
for file in Path("output").glob("*_final.parquet"):
    gene = file.stem.split("_")[0]
    dfs[gene] = pl.read_parquet(f"output/{gene}_final.parquet")
    sams[gene] = pl.read_parquet(f"output/{gene}_all.parquet")
    filtereds[gene] = pl.read_parquet(f"output/{gene}_filtered.parquet")
    offtargets[gene] = pl.read_parquet(f"output/{gene}_offtargets.parquet")
    try:
        overlapped[gene] = pl.read_parquet(f"output/{gene}_final_overlapped.parquet")
    except FileNotFoundError:
        pass
# %%

sams["Abi3"].filter(
    pl.col("match_max").is_between(pl.col("length") * 0.5, pl.col("length") * 0.8)
    & pl.col("match_max").gt(15)
).with_columns(
    tm_match=pl.struct(["seq", "cigar", "mismatched_reference"]).apply(
        lambda x: tm_match(x["seq"], x["cigar"], x["mismatched_reference"])
    )
).filter(
    pl.col("tm_match").gt(40)
)

# %%
for x in (
    sams["Abi3"]
    .filter(
        pl.col("match_max").is_between(pl.col("length") * 0.6, pl.col("length") * 0.7)
        & pl.col("match_max").gt(16)
    )
    .iter_rows(named=True)
):
    print(x["id"])
    print(tm_match(x["seq"], x["cigar"], x["mismatched_reference"]))
# %%
length = {k: len(v) for k, v in dfs.items()}
# filter length less than 30
short = {k: v for k, v in length.items() if v < 30}
print(len(length), len(short))
# %%
stripplot(all=sam["pos_end"], filtered=filtered["pos_end"], s=df["pos_end"])


# %%
def runpls(gene: str):
    subprocess.run(
        ["python", "scripts/new_postprocess.py", gene, "-O", "10"],
        check=True,
        capture_output=True,
    )


with ThreadPoolExecutor(32) as executor:
    for x in as_completed([executor.submit(runpls, gene) for gene in short.keys()]):
        print("ok")
        x.result()


# %%

# %%

model = gen_model(47)
g = secondary_structure(
    "TTGCTAGCGTGGGCCAATTAGAGTGAGTAGTAGTGGAGTAAGTAGCTGGGACCCTGAAGCTGGTCCCATAAAGAGTGAGTAGTAGTGGAGTATGTGGTTTGGAGATGATAGACCCTATAGTGAGTCGTATTAGTTGTCAGCCATTGCGGG"
)

# %%
g = nonspecific_test(
    "TAGTGGAGTAAGTAGCTGGGACCCTGAAGC",
    "TTGCTAGCGTGGGCCAATTAGAGTGAGTAGTAGTGGAGTAAGTAGCTGGGACCCTGAAGCTGGTCCCATAAAGAGTGAGTAGTAGTGGAGTATGTGGTTTGGAGATGATAGACCCTATAGTGAGTCGTATTAGTTGTCAGCCATTGCGGG",
    t=47,
)

# %%
tm_fish
sam.with_columns(mapped_to=pl.col("transcript").apply(gtf_all.ts_to_gene))
# %%
probe = reverse_complement("ATTTCAGTGGCTGCTACTCGGCGCTT")
seq = gtf_all.get_seq("ENSMUST00000134932")

# %%
nonspecific_test(probe, seq[1147:1177], t=47)

# %%
gtf_all.get_seq("ENSMUST00000102476")
# %%


# %%
# %%
