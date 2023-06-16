# %%
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import cycle, permutations
from pathlib import Path
from typing import Iterable, Iterator, Sequence

import numpy as np
import polars as pl

from oligocheck.geneframe import GeneFrame
from oligocheck.merfish.alignment import gen_fasta, gen_fastq, run_bowtie
from oligocheck.external.external_data import ExternalData
from oligocheck.merfish.pairwise import pairwise_alignment
from oligocheck.picker.codebook import CodebookPickerSingleCell
from oligocheck.seqcalc import tm_hybrid
from oligocheck.sequtils import reverse_complement, slide

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

# ori = pl.read_excel("panels/motorcortex_ori.xlsx", read_csv_options={"skip_rows": 1})
# orir = pl.read_excel("panels/motorcortex_ori_readout.xlsx", read_csv_options={"skip_rows": 1})
# sm = orir.filter(~pl.col("Purpose").str.starts_with("MERFISH"))["Purpose"]
# # %%
# # %%
# convert = json.loads(Path("panels/motorcortex_convert.json").read_text())
# convert["01-Mar"] = "Marchf1"

# mapping = dict(orir[["Sequence", "Readout probe name"]].iter_rows())
# codebook = (
#     ori.filter(~pl.col("Gene name").is_in(sm))
#     .with_columns(
#         ro1=pl.col("Sequence").str.split(" ").list.get(2),
#         ro2=pl.col("Sequence").str.split(" ").list.get(5),
#     )
#     .select(["Gene name", "ro1", "ro2"])
#     .groupby("Gene name")
#     .all()
#     .select("Gene name", combi=pl.col("ro1").list.concat(pl.col("ro2")).list.unique())
#     .explode("combi")
#     .with_columns(
#         pl.col("Gene name").apply(lambda x: convert.get(x, x)).alias("Gene name"),
#         pl.col("combi").apply(reverse_complement).apply(mapping.get).alias("combi"),
#     )
# )
# bitnmap = dict(
#     codebook[["combi"]]
#     .unique("combi")
#     .sort("combi")
#     .with_row_count("id", offset=1)[["combi", "id"]]
#     .iter_rows()
# )

# codebook = dict(
#     codebook.with_columns(combi=pl.col("combi").apply(bitnmap.get)).groupby("Gene name").all().iter_rows()
# )

# blacklist = set(pl.read_csv("data/readout_fused_bad.csv")[["split1", "split2"]].iter_rows())
# codebook = pl.read_csv("panels/celltype_codebook.csv")
name = 'sm'
readouts = pl.read_csv("data/readout_ref_filtered.csv")
# genes = Path("zach_26.txt").read_text().splitlines()
# genes = Path(f"{name}_25.txt").read_text().splitlines()
smgenes = ['Cdc42','Neurog2','Ccnd2']
genes, _, _ = gtf_all.check_gene_names(smgenes)
acceptable_tss = {g: set(pl.read_csv(f"output/{g}_acceptable_tss.csv")["transcript"]) for g in genes}
n, short_threshold = 67, 65
# %%
dfx, overlapped = {}, {}
for gene in genes:
    dfx[gene] = GeneFrame.read_parquet(f"output/{gene}_final.parquet")
dfs = GeneFrame.concat(dfx.values())
short = dfs.count("gene").filter(pl.col("count") < short_threshold)


# %%
fixed_n = {}
short_fixed = {}


def run_overlap(genes: Iterable[str], overlap: int):
    def runpls(gene: str):
        subprocess.run(
            ["python", "scripts/new_postprocess.py", gene, "-O", str(overlap)],
            check=True,
            capture_output=True,
        )

    with ThreadPoolExecutor(32) as executor:
        for x in as_completed(
            [
                executor.submit(runpls, gene)
                for gene in genes
                if not Path(f"output/{gene}_final_overlap_{overlap}.parquet").exists()
            ]
        ):
            print("ok")
            x.result()


needs_fixing = set(short["gene"])

for ol in [5, 10, 15, 20]:
    print(ol, needs_fixing)
    run_overlap(needs_fixing, ol)
    for gene in needs_fixing.copy():
        df = GeneFrame.read_parquet(f"output/{gene}_final_overlap_{ol}.parquet")
        if len(df) >= short_threshold or ol == 20:
            needs_fixing.remove(gene)
            fixed_n[gene] = ol
            short_fixed[gene] = df
# else:
#     raise ValueError(f"Gene {gene} cannot be fixed")

short_fixed = GeneFrame.concat(short_fixed.values())
# %%
cutted = GeneFrame.concat([dfs.filter(~pl.col("gene").is_in(short["gene"])), short_fixed[dfs.columns]])
cutted = GeneFrame(
    cutted.sort(["gene", "priority"]).groupby("gene").agg(pl.all().head(n)).explode(pl.all().exclude("gene"))
)
# %%
counts = cutted.count("gene")
print(len(cutted))

print("dumped", counts.filter(pl.col("count") < 40))
counts = counts.filter(pl.col("count") >= 40)

# %%
scdata = pl.read_parquet("data/fpkm/embryo_raw.parquet")

# %%
# cl = counts[:, perc999.filter(pl.col("counts") < 50)["gene"]].to_numpy()
import importlib

import oligocheck.picker.codebook as CB

importlib.reload(CB)

if name == "zach":

    code_existing = json.loads(Path("constructed/cs.json").read_text())
    code_matrix = np.zeros((len(code_existing), 27), dtype=bool)
    for i, row in enumerate(code_existing.values()):
        for code in row:
            code_matrix[i, code-1] = 1
    assert code_matrix[:,24:].sum() == 0
    code_gene = scdata[:, list(code_existing.keys())]
    counts = counts.filter(~pl.col("gene").is_in(list(code_existing.keys())))

    mhd = CB.CodebookPickerSingleCell("static/mhd4_27bit.csv", code_existing=code_matrix, counts_existing=code_gene.to_numpy())
    best, perro = mhd.find_optimalish(scdata[:, counts["gene"]], iterations=100)
    codebook = mhd.gen_codebook(int(best))[: counts.shape[0]]
    assert not set(tuple(x) for x in codebook) & set(tuple(x) for x in code_matrix)
    codebook = np.where(codebook)[1].reshape([-1, 4])
    codebook = dict(zip(counts['gene'], (codebook+1).tolist()))
elif name =="cs":
    mhd = CB.CodebookPickerSingleCell("static/mhd4_27bit.csv", subset=(0,24))
    best, perro = mhd.find_optimalish(scdata[:, counts["gene"]], iterations=100)
    codebook = mhd.gen_codebook(int(best))[: counts.shape[0]]
    codebook = np.where(codebook)[1].reshape([-1, 4])
    codebook = dict(zip(counts['gene'], (codebook+1).tolist()))
elif name == "sm":
    codebook = {g:[i] for i,g in enumerate(genes, 28)}


#%%
# assert not codebook[:,24:].sum()




# %%
def stitch(seq: str, codes: Sequence[str], sep: str = "TT") -> str:
    return codes[0] + "TT" + codes[0] + sep + seq + sep + codes[1]


def gen_stitch(n: int, seq: str, seq_map: dict[str, str], ros: Iterator[tuple[str, str]]):
    for _ in range(n):
        codes = next(ros)
        for sep in ["TT", "AA", "AT", "TA"]:
            consd = stitch(seq, [seq_map[x] for x in codes], sep=sep)
            if not (set(slide(consd, 15)) & trna_rna_kmers) and not (set(slide(consd, 18)) & kmerset):
                yield consd, codes


def construct_encoding(seq_encoding: pl.DataFrame, return_all: bool = False, sm:bool=False):
    gene = seq_encoding[0, "gene"]
    ros = readouts.filter(pl.col("id").is_in(codebook[gene]))
    if sm and not len(codebook[gene]) == 1:
        raise ValueError(f"Gene {gene} has multiple readouts but sm is True")
    seq_map = dict(ros[["name", "seq"]].iter_rows())

    # No blacklists here.
    if not sm:
        fusedreadouts: list[tuple[str, str]] = [cs for cs in permutations(ros["name"], 2)]
    else:
        fusedreadouts: list[tuple[str, str]] = [(ros[0, 'name'], ros[0, 'name'])]

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

            out["name"].append(name + "" if not return_all else f"{name}__{i}")
            out["constructed"].append(reverse_complement(consd))
            out["code1"].append(codes[0])
            out["code2"].append(codes[1])
            if not return_all:
                break

    return pl.DataFrame(out)


headers = {"cs": "TGGCGACTGAAGTACGAGTCC", "zach": "GAGAGGCGAGGACACCTACAG", 'sm': "CAACCGTACCGCTTGCTTACC"}


# %%
constructed = []
for gene, df_ in cutted.filter(pl.col("gene").is_in(counts['gene'])).groupby("gene"):
    df_ = df_.join(construct_encoding(df_, sm=name == "sm"), on="name", how="inner").with_columns(constructed=headers[name] + pl.col("constructed") + "TATTTCCC")
    constructed.append(df_)
constructed = GeneFrame.concat(constructed)
# constructed = (
#     cutted.groupby("gene")
#     .apply(lambda df: df.join(construct_encoding(df), on="name", how="inner"))
#     .with_columns(constructed=headers["celltype"] + pl.col("constructed") + "TATTTCCC")
# )


# %%


def check_offtargets(constructed: pl.DataFrame, acceptable_tss: dict[str, list[str]]):
    out = []
    nogos = []
    for gene, df in GeneFrame.from_sam(
        run_bowtie(
            gen_fastq(constructed["name"], constructed["constructed"]).getvalue(),
            "data/mm39/mm39",
            seed_length=15,
            threshold=20,
            n_return=500,
        )
    ).groupby("gene"):
        filtered, nogo = GeneFrame(df).filter_by_match(acceptable_tss[gene], match_max=0.8, return_nogo=True)
        out.append(filtered)
        nogos += nogo

    return nogos
    # tm_offtarget = (
    #     y.filter(
    #         pl.col("match_max").is_between(pl.col("length") * 0.5, pl.col("length") * 0.8 + 0.01)
    #         & pl.col("match_max").gt(15)
    #     )
    #     .with_columns(
    #         tm_offtarget=pl.struct(["seq", "cigar", "mismatched_reference"]).apply(
    #             lambda x: tm_match(x["seq"], x["cigar"], x["mismatched_reference"])
    #         )
    #     )
    #     .groupby("name")
    #     .agg(pl.max("tm_offtarget"))
    #     .filter(pl.col("tm_offtarget").gt(40))
    # )
    # offtarget = (
    #     y.filter(
    #         pl.col("match_max").gt(0.8 * pl.col("length"))
    #         # & pl.col("gene").ne(pl.col("transcript_name"))
    #         & pl.col("flag").map(lambda x: x & 16 > 0)
    #     )
    #     .groupby("gene")
    #     .apply(lambda df: df.filter(~pl.col("transcript").is_in(acceptable_tss[df["gene"][0]])))
    #     .groupby("name")
    #     .agg(pl.max("match_max"))
    # )

    # return offtarget["name"].to_list() + tm_offtarget["name"].to_list()


nogo = check_offtargets(constructed, acceptable_tss)
pass1 = constructed.filter(~pl.col("name").is_in(nogo))
len(nogo)
# %%

pass2 = []

for _, df in cutted.filter(pl.col("name").is_in(nogo)).groupby("gene"):
    pass2.append(
        construct_encoding(df, sm=name=="sm", return_all=True).with_columns(
            constructed=headers[name] + pl.col("constructed") + "TATTTCCC"
        )
    )
pass2 = pl.concat(pass2).with_columns(
    name_var=pl.col("name"),
    name=pl.col("name").str.split("__").list.get(0),
)
# %%
nogo2 = check_offtargets(pass2, acceptable_tss)
# %%
constructed = GeneFrame.concat(
    [
        pass1,
        cutted.join(
            pass2.filter(~pl.col("name_var").is_in(nogo2)).unique("name").drop("name_var"),
            on="name",
            how="inner",
        ),
    ]
)

len(constructed)
# %%

footers = {"zach": "TATAGTGAGTCGTATTAGACCGGTCT", "cs": "TATAGTGAGTCGTATTAGAGGCACTG", 'sm': 'TATAGTGAGTCGTATTAAGGCGGTT' }
constructed.with_columns(constructed=pl.col("constructed") + footers[name]).write_parquet(f"constructed/{name}.parquet")

Path(f"constructed/{name}.json").write_text(json.dumps(codebook, indent=2))
#%%
#%%
combi =GeneFrame.concat([pl.read_parquet("constructed/zach.parquet"), pl.read_parquet("constructed/cs.parquet"), pl.read_parquet("constructed/sm.parquet")])

len(combi)

# %%
combi.write_parquet("constructed/combi.parquet")
combi[['gene']].unique('gene').sort('gene').write_csv("constructed/combi.txt", has_header=False)
#%%


def backfill(seq: str, target:int= 167):
    return "GGTAAGAAGTGAGTTAGTGGAA"[:max(0, target - len(seq))] + seq

combi_filled = combi.with_columns(constructed=pl.col("constructed").apply(backfill))

t7 = reverse_complement('TAATACGACTCACTATAGGGAAATA')
assert combi_filled.select(pl.col("constructed").str.contains(t7))['constructed'].all()
assert combi_filled.select(pl.col("constructed").str.contains("|".join(footers.values())))['constructed'].all()
assert combi_filled.select(pl.col("constructed").str.contains("|".join(headers.values())))['constructed'].all()

length = combi_filled['constructed'].apply(len)
print(length.min(), length.max())
combi_filled[['constructed']].write_csv("constructed/combi_filled.txt", has_header=False)

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

theirs = (
    ori.filter(~pl.col("Gene name").is_in(sm))
    .with_row_count("id")
    .with_columns(seq=pl.col("Sequence").str.split(" ").list.get(3))
    .with_columns(name=pl.col("Gene name") + "_" + pl.col("id").cast(pl.Utf8))
)

res = GeneFrame.from_sam(
    run_bowtie(
        gen_fasta(
            theirs["name"], theirs.select(pl.col("Sequence").str.replace(" ", "", n=10))["Sequence"]
        ).getvalue(),
        "data/mm39/mm39",
        fasta=True,
    ),
    split_name=False,
).with_columns(gene=pl.col("name").str.split("_").list.get(0))
# %%
ontarget = res.filter_eq(aln_score=60).unique("name").with_columns(tm=pl.col("seq").apply(tm_hybrid))
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
))