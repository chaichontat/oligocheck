# %%
from itertools import cycle, permutations
from pathlib import Path
from typing import Sequence

import numpy as np
import polars as pl

from oligocheck.merfish.alignment import gen_fastq, run_bowtie
from oligocheck.merfish.external_data import ExternalData
from oligocheck.merfish.filtration import count_match
from oligocheck.merfish.pairwise import pairwise_alignment
from oligocheck.sequtils import parse_sam, reverse_complement, stripplot

gtf_all = ExternalData(
    cache="data/mm39/gencode_vM32_transcripts_all.parquet",
    path="data/mm39/Mus_musculus.GRCm39.109.gtf",
    fasta="data/mm39/combi.fa.gz",
)


readouts = pl.read_csv("data/readout_ref_filtered.csv")
codebook = pl.read_csv("panels/celltype_codebook.csv")
blacklist = set(pl.read_csv("data/readout_fused_bad.csv")[["split1", "split2"]].iter_rows())
genes = Path("celltype_genes.csv").read_text().splitlines()
acceptable_tss = {g: set(pl.read_csv(f"output/{g}_acceptable_tss.csv")["transcript"]) for g in genes}
n = 40
# %%
dfs, overlapped = {}, {}
for gene in genes:
    dfs[gene] = pl.read_parquet(f"output/{gene}_final.parquet")
    try:
        overlapped[gene] = pl.read_parquet(f"output/{gene}_final_overlapped.parquet")
    except FileNotFoundError:
        pass


# %%

dfs = pl.concat((dfs | overlapped).values())
# %%
counts = dfs.groupby("transcript_name").agg(pl.count("pos_end"))


# %%
def stitch(seq: str, codes: Sequence[str]) -> str:
    return codes[0] + "TT" + codes[1] + "TT" + seq + "TT" + codes[2] + "TT" + codes[3]


def construct_encoding(seq_encoding: pl.DataFrame):
    ros = readouts.filter(pl.col("id").is_in(np.nonzero(codebook[0, 1:])[1] + 1))
    seq_map = dict(ros[["name", "seq"]].iter_rows())

    # Filter out blacklists
    fusedreadouts: list[tuple[str, str, str, str]] = [
        cs
        for cs in permutations(ros["name"], 4)
        if not ((cs[0], cs[1]) in blacklist or (cs[2], cs[3]) in blacklist)
    ]
    assert len(fusedreadouts) > 4

    out = dict(name=[], constructed=[], code1=[], code2=[], code3=[], code4=[])

    # rand.shuffle(combi)
    for (name, seq), codes in zip(seq_encoding[["name", "seq"]].iter_rows(), cycle(fusedreadouts)):
        # Reverse complemented here.
        # Always use mRNA sequence and readouts with Cs.
        out["name"].append(name)
        out["constructed"].append(reverse_complement(stitch(seq, [seq_map[x] for x in codes])))
        out["code1"].append(codes[0])
        out["code2"].append(codes[1])
        out["code3"].append(codes[2])
        out["code4"].append(codes[3])

    return seq_encoding.join(pl.DataFrame(out), on="name", how="left")


headers = {"celltype": "GAGAGGCGAGGACACCTACAG"}
footers = {"celltype": "TATTTCCCTATAGTGAGTCGTATTAGACCGGTCT"}

# %%
constructed = (
    dfs.groupby("transcript_name")
    .apply(construct_encoding)
    .with_columns(constructed=headers["celltype"] + pl.col("constructed") + "TATTTCCC")
)

# %%
y = count_match(
    parse_sam(
        run_bowtie(
            gen_fastq(constructed["name"], constructed["constructed"]).getvalue(),
            "data/mm39/mm39",
            seed_length=15,
            threshold=20,
            n_return=500,
        )
    )
)

y
# %%
offtarget = (
    y.filter(
        pl.col("match_max").gt(0.75 * pl.col("length"))
        # & pl.col("gene").ne(pl.col("transcript_name"))
        & pl.col("flag").map(lambda x: x & 16 > 0)
    )
    .groupby("gene")
    .apply(lambda df: df.filter(~pl.col("transcript").is_in(acceptable_tss[df["gene"][0]])))
    .groupby("name")
    .agg(pl.max("match_max"))
)

len(constructed), len(offtarget)

# %%
finalfiltered = constructed.filter(~pl.col("name").is_in(offtarget["name"]))

# .groupby("gene").agg(pl.count()).sort("count")

# %%
constructed.groupby("gene").agg(pl.count()).sort("count").filter(pl.col("count").lt(20))


# %%
def extract_seqs(df: pl.DataFrame, n: int | None = 40):
    return df.sort("maps_to_pseudo", "priority", "pos_end", descending=[False, False, True])[:n]


# seqs = {k: extract_seqs(v, n) for k, v in dfs.items()}
# seqs |= {k: extract_seqs(v, None) for k, v in overlapped.items()}
# filtered = pl.concat(list(seqs.values()))
# %%


p = pairwise_alignment(y[5, "seq"], y[5, "cigar"], y[5, "mismatched_reference"])
[tm_hybrid(y[5, "seq"][slice(*x.span())]) for x in r.finditer(p[1])]

# %%
