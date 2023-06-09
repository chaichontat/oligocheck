# %%
import json
from itertools import cycle, permutations
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import polars as pl

from oligocheck.geneframe import GeneFrame
from oligocheck.merfish.alignment import gen_fasta, gen_fastq, run_bowtie
from oligocheck.merfish.external_data import ExternalData
from oligocheck.merfish.filtration import count_match
from oligocheck.merfish.pairwise import pairwise_alignment
from oligocheck.seqcalc import tm_hybrid, tm_match
from oligocheck.sequtils import parse_sam, reverse_complement, slide, stripplot

gtf_all = ExternalData(
    cache="data/mm39/gencode_vM32_transcripts_all.parquet",
    path="data/mm39/Mus_musculus.GRCm39.109.gtf",
    fasta="data/mm39/combi.fa.gz",
)
kmer18 = pl.read_csv(
    "data/mm39/kmer_genome18.txt", separator=" ", has_header=False, new_columns=["kmer", "count"]
)
# rrnas = get_rrna("data/mm39/Mus_musculus.GRCm39.ncrna.fa.gz")
trna_rna_kmers = set(
    pl.read_csv(
        "data/mm39/kmer_trcombi15.txt", separator=" ", has_header=False, new_columns=["kmer", "count"]
    )["kmer"]
)


kmerset = set(kmer18["kmer"])

ori = pl.read_excel("panels/motorcortex_ori.xlsx", read_csv_options={"skip_rows": 1})
orir = pl.read_excel("panels/motorcortex_ori_readout.xlsx", read_csv_options={"skip_rows": 1})
sm = orir.filter(~pl.col("Purpose").str.starts_with("MERFISH"))["Purpose"]
# %%

theirs = (
    ori.filter(~pl.col("Gene name").is_in(sm))
    .with_row_count("id")
    .with_columns(seq=pl.col("Sequence").str.split(" ").arr.get(3))
    .with_columns(name=pl.col("Gene name") + "_" + pl.col("id").cast(pl.Utf8))
)

res = count_match(
    parse_sam(
        run_bowtie(
            gen_fasta(
                theirs["name"], theirs.select(pl.col("Sequence").str.replace(" ", "", n=10))["Sequence"]
            ).getvalue(),
            "data/mm39/mm39",
            fasta=True,
        ),
        split_name=False,
    )
).with_columns(gene=pl.col("name").str.split("_").arr.get(0))
# %%
ontarget = (
    res.filter(pl.col("aln_score") == 60).unique("name").with_columns(tm=pl.col("seq").apply(tm_hybrid))
)
offtarget = (
    res.filter(pl.col("match_max") > 24)
    .with_columns(transcript_name=pl.col("transcript").apply(gtf_all.ts_to_gene))
    .filter(pl.col("gene") != pl.col("transcript_name"))
    .filter(~pl.col("transcript_name").str.starts_with("Gm"))
    .sort("match_max")
    .unique("name", keep="first")
)
# %%
offtarget.groupby("gene").agg(pl.max("match_max"))
b = (
    offtarget.filter(pl.col("match_max") > 24)
    .with_columns(transcript_name=pl.col("transcript").apply(gtf_all.ts_to_gene))
    .filter(pl.col("gene") != pl.col("transcript_name"))
    .filter(~pl.col("transcript_name").str.starts_with("Gm"))
)
# %%
convert = json.loads(Path("panels/motorcortex_convert.json").read_text())
convert["01-Mar"] = "Marchf1"

mapping = dict(orir[["Sequence", "Readout probe name"]].iter_rows())
codebook = (
    ori.filter(~pl.col("Gene name").is_in(sm))
    .with_columns(
        ro1=pl.col("Sequence").str.split(" ").arr.get(2),
        ro2=pl.col("Sequence").str.split(" ").arr.get(5),
    )
    .select(["Gene name", "ro1", "ro2"])
    .groupby("Gene name")
    .all()
    .select("Gene name", combi=pl.col("ro1").arr.concat(pl.col("ro2")).arr.unique())
    .explode("combi")
    .with_columns(
        pl.col("Gene name").apply(lambda x: convert.get(x, x)).alias("Gene name"),
        pl.col("combi").apply(reverse_complement).apply(mapping.get).alias("combi"),
    )
)
bitnmap = dict(
    codebook[["combi"]]
    .unique("combi")
    .sort("combi")
    .with_row_count("id", offset=1)[["combi", "id"]]
    .iter_rows()
)

# %%
readouts = pl.read_csv("data/readout_ref_filtered.csv")
codebook = dict(
    codebook.with_columns(combi=pl.col("combi").apply(bitnmap.get)).groupby("Gene name").all().iter_rows()
)

# codebook = pl.read_csv("panels/celltype_codebook.csv")
# blacklist = set(pl.read_csv("data/readout_fused_bad.csv")[["split1", "split2"]].iter_rows())
genes = Path("panels/motorcortex_converted.txt").read_text().splitlines()
acceptable_tss = {g: set(pl.read_csv(f"output/{g}_acceptable_tss.csv")["transcript"]) for g in genes}
n = 48
# %%
dfs, overlapped = {}, {}
for gene in genes:
    dfs[gene] = pl.read_parquet(f"output/{gene}_final.parquet")
dfs = pl.concat(dfs.values())
short = dfs.groupby("gene").agg(pl.count()).filter(pl.col("count") < 45)


# %%
fixed_n = {}
short_fixed = {}
for gene in short["gene"]:
    for ol in [5, 10, 15, 20]:
        df = pl.read_parquet(f"output/{gene}_final_overlap_{ol}.parquet")
        fixed_n[gene] = ol
        if len(df) >= n:
            short_fixed[gene] = df
            break
        if ol == 20:
            short_fixed[gene] = df
            break
    else:
        raise ValueError(f"Gene {gene} cannot be fixed")
short_fixed = pl.concat(short_fixed.values())

cutted = (
    GeneFrame(dfs)
    .sort(["gene", "priority"])
    .groupby("gene")
    .agg(pl.all().head(48))
    .explode(pl.all().exclude("gene"))
)
dfs = pl.concat([cutted.filter(~pl.col("gene").is_in(short["gene"])), short_fixed[cutted.columns]])
# %%
counts = dfs.groupby("gene").agg(pl.count("pos_end"))


# %%
def stitch(seq: str, codes: Sequence[str], sep: str = "TT") -> str:
    return codes[0] + "TT" + codes[0] + sep + seq + sep + codes[1] + "TT" + codes[1]


def gen_stitch(n, seq: str, seq_map: dict[str, str], ros: Iterable[tuple[str, str]]):
    for _ in range(n):
        codes = next(ros)
        for sep in ["TT", "AA", "AT", "TA"]:
            consd = stitch(seq, [seq_map[x] for x in codes], sep=sep)
            if not (set(slide(consd, 15)) & trna_rna_kmers) and not (set(slide(consd, 18)) & kmerset):
                yield consd, codes


def construct_encoding(seq_encoding: pl.DataFrame, return_all: bool = False):
    gene = seq_encoding[0, "gene"]
    ros = readouts.filter(pl.col("id").is_in(codebook[gene]))
    seq_map = dict(ros[["name", "seq"]].iter_rows())

    # No blacklists here.
    fusedreadouts: list[tuple[str, str]] = [cs for cs in permutations(ros["name"], 2)]

    out = dict(name=[], constructed=[], code1=[], code2=[])
    ros = cycle(fusedreadouts)
    # rand.shuffle(combi)
    for name, seq in seq_encoding[["name", "seq"]].iter_rows():
        # Reverse complemented here.
        # Always use mRNA sequence and readouts with Cs.
        generator = gen_stitch(len(fusedreadouts), seq, seq_map, ros)
        i = 0
        while True:
            try:
                consd, codes = next(generator)
                i += 1
            except StopIteration:
                # problematic = (set(slide(consd, 18)) & kmerset) | (set(slide(consd, 15)) & trna_rna_kmers)
                # print(problematic, [consd.find(x) for x in problematic])
                # print(f"Cannot find a proper readout for {name}")
                break
            # print("consd")

            out["name"].append(name + "" if not return_all else f"{name}__{i}")
            out["constructed"].append(reverse_complement(consd))
            out["code1"].append(codes[0])
            out["code2"].append(codes[1])
            if not return_all:
                break

    return pl.DataFrame(out)


headers = {"celltype": "GAGAGGCGAGGACACCTACAG"}
footers = {"celltype": "TATTTCCCTATAGTGAGTCGTATTAGACCGGTCT"}

# %%
constructed = (
    dfs.groupby("gene")
    .apply(lambda df: df.join(construct_encoding(df), on="name", how="inner"))
    .with_columns(constructed=headers["celltype"] + pl.col("constructed") + "TATTTCCC")
    # .filter(~pl.col("constructed").apply(lambda x: bool(set(slide(x, 18)) & kmerset)))
    # .filter(~pl.col("constructed").apply(lambda x: bool(set(slide(x, 15)) & trna_rna_kmers)))
)


# %%


def check_offtargets(constructed: pl.DataFrame, acceptable_tss: dict[str, list[int]]):
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
    tm_offtarget = (
        y.filter(
            pl.col("match_max").is_between(pl.col("length") * 0.5, pl.col("length") * 0.8 + 0.01)
            & pl.col("match_max").gt(15)
        )
        .with_columns(
            tm_offtarget=pl.struct(["seq", "cigar", "mismatched_reference"]).apply(
                lambda x: tm_match(x["seq"], x["cigar"], x["mismatched_reference"])
            )
        )
        .groupby("name")
        .agg(pl.max("tm_offtarget"))
        .filter(pl.col("tm_offtarget").gt(40))
    )
    offtarget = (
        y.filter(
            pl.col("match_max").gt(0.8 * pl.col("length"))
            # & pl.col("gene").ne(pl.col("transcript_name"))
            & pl.col("flag").map(lambda x: x & 16 > 0)
        )
        .groupby("gene")
        .apply(lambda df: df.filter(~pl.col("transcript").is_in(acceptable_tss[df["gene"][0]])))
        .groupby("name")
        .agg(pl.max("match_max"))
    )

    return offtarget["name"].to_list() + tm_offtarget["name"].to_list()


nogo = check_offtargets(constructed, acceptable_tss)
pass1 = constructed.filter(~pl.col("name").is_in(nogo))
# %%

pass2 = []

for _, df in dfs.filter(pl.col("name").is_in(nogo)).groupby("gene"):
    pass2.append(
        construct_encoding(df, return_all=True).with_columns(
            constructed=headers["celltype"] + pl.col("constructed") + "TATTTCCC"
        )
    )
pass2 = pl.concat(pass2).with_columns(
    name_var=pl.col("name"),
    name=pl.col("name").str.split("__").arr.get(0),
)
# %%
nogo2 = check_offtargets(pass2, acceptable_tss)
# %%
constructed = pl.concat(
    [
        pass1,
        dfs.join(
            pass2.filter(~pl.col("name_var").is_in(nogo2)).unique("name").drop("name_var"),
            on="name",
            how="inner",
        ),
    ]
)

len(constructed)
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
# %%
