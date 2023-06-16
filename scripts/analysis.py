# %%
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import polars as pl

from oligocheck.external.external_data import ExternalData
from oligocheck.geneframe import GeneFrame
from oligocheck.merfish.alignment import gen_fasta
from oligocheck.merfish.nnupack import gen_model, nonspecific_test, secondary_structure
from oligocheck.seqcalc import tm_fish, tm_match
from oligocheck.sequtils import reverse_complement, stripplot

gtf = ExternalData(
    cache="data/mm39/gencode_vM32_transcripts.parquet",
    path="data/mm39/gencode.vM32.chr_patch_hapl_scaff.basic.annotation.gtf",
    fasta="data/mm39/combi.fa.gz",
)


gtf_all = ExternalData(
    cache="data/mm39/gencode_vM32_transcripts_all.parquet",
    path="data/mm39/Mus_musculus.GRCm39.109.gtf",
    fasta="data/mm39/combi.fa.gz",
)
pl.Config.set_fmt_str_lengths(100)
pl.Config.set_tbl_rows(100)
# %%
snos = gtf_all.gtf.filter(pl.col("gene_biotype") == "snoRNA")["transcript_id"]
df = pl.DataFrame(dict(name=snos, seq=snos.apply(gtf_all.get_seq)))
Path("data/snorna.fa").write_text(gen_fasta(names=df["name"], seqs=df["seq"]).getvalue())
# %%
genes = Path("panels/motorcortex_smgenes.txt").read_text().splitlines()
# %%
dfx, sams, filtereds, offtargets, overlapped = {}, {}, {}, {}, {}
for gene in genes:
    dfx[gene] = pl.read_parquet(f"output/{gene}_final.parquet")
    sams[gene] = pl.read_parquet(f"output/{gene}_all.parquet")
    filtereds[gene] = pl.read_parquet(f"output/{gene}_filtered.parquet")
    offtargets[gene] = pl.read_parquet(f"output/{gene}_offtargets.parquet")
    # try:
    #     overlapped[gene] = pl.read_parquet(f"output/{gene}_final_overlapped.parquet")
    # except FileNotFoundError:
    #     pass
dfx = dfx | overlapped
# %%
dfs = GeneFrame.concat(dfx.values()).sort(
    ["gene", pl.col("isoforms").list.lengths(), "priority", "pos_end"],
    descending=[False, True, False, True],
)

# cutted = (
#     GeneFrame(pl.concat(dfs.values()).sort(["gene", "priority"]))
#     .groupby("gene")
#     .agg(pl.all().head(48))
#     .explode(pl.all().exclude("gene"))
# )

# %%
counts = dfs.count()
# length = {k: len(v) for k, v in counts.iter_rows()}
# filter length less than 30
short = {k: v for k, v in counts.iter_rows() if v < 48}
print(len(counts), len(short))
# %%
target = "Mup5"
stripplot(
    all=sams[target]["pos_end"],
    filtered=filtereds[target]["pos_end"],
    selected=dfs.gene(target)["pos_end"],
)


# %%
# %%


# %%
def runpls(gene: str):
    subprocess.run(
        ["python", "scripts/new_postprocess.py", gene, "-O", "20"],
        check=True,
        capture_output=True,
    )


with ThreadPoolExecutor(32) as executor:
    for x in as_completed([executor.submit(runpls, gene) for gene in short.keys()]):
        print("ok")
        x.result()


# %%
fixed_n = {}
short_fixed = {}
for gene in short.keys():
    for ol in [5, 10, 15, 20]:
        df = pl.read_parquet(f"output/{gene}_final_overlap_{ol}.parquet")
        fixed_n[gene] = ol
        if len(df) >= 48:
            short_fixed[gene] = df
            break
        if ol == 20:
            short_fixed[gene] = df

    # shortfix.groupby("gene").agg(pl.count())
# short_fixed = pl.concat(short_fixed.values())

# dfs = pl.concat([cutted.filter(~pl.col("gene").is_in(short.keys()))[short_fixed.columns], short_fixed])
# %%
stripplot(all=sam["pos_end"], filtered=filtered["pos_end"], s=df["pos_end"])

# %%
from oligocheck.merfish.nnupack import gen_model, nonspecific_test, secondary_structure

dfs.lazy().gene("Slc17a6").with_columns(
    secondary=pl.col("seq").apply(lambda x: secondary_structure(x)["(probe)"][3])
)
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
