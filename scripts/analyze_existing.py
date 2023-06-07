# %%
import polars as pl

d = (
    pl.read_csv("data/published_celltype.csv")
    .with_row_count("id")
    .with_columns(
        name=pl.col("Gene") + "_" + pl.col("id").cast(pl.Utf8),
        binding=pl.col("Probe Sequence").str.slice(40, 30),
    )
    .with_columns(
        tm=pl.col("binding").apply(tm_fish) - 0.65 * 30,
        hp=pl.col("binding").apply(hp_fish) - 0.65 * 30,
    )
)
# %%
res = parse_sam(
    run_bowtie(gen_fastq(d.with_columns(seq=pl.col("binding"))).getvalue(), "data/mm39/mm39"),
    split_name=False,
).with_columns(
    transcript_name=pl.col("transcript").apply(gtf_all.ts_to_gene),
)
# %%
t = res.filter(pl.col("name").str.starts_with("Abi3") & pl.col("aln_score").eq(60)).sort("aln_score")

# %%
res.filter(pl.col("name") == "Abi3_88").with_columns(
    transcript=pl.col("transcript").str.extract(r"(ENSMUST\d+)")
).join(fpkm, left_on="transcript", right_on="transcript_id(s)", how="left")
# %%
from Levenshtein import distance
import polars as pl

ro = pl.read_csv("data/readout_ref.csv")

m = 1000
for i in range(len(ro)):
    for j in range(i + 1, len(ro)):
        m = min(distance(ro[i, "Sequence"], ro[j, "Sequence"]), m)
print(m)

def gen
# %%

d = (
    pl.read_csv("data/published_celltype.csv")
    .with_row_count("id")
    .with_columns(
        name=pl.col("Gene") + "_" + pl.col("id").cast(pl.Utf8),
        binding=pl.col("Probe Sequence").str.slice(40, 30),
    )
    .with_columns(
        tm=pl.col("binding").apply(tm_fish) - 0.65 * 30,
        hp=pl.col("binding").apply(hp_fish) - 0.65 * 30,
    )
)
# %%
res = parse_sam(
    run_bowtie(gen_fastq(d.with_columns(seq=pl.col("binding"))).getvalue(), "data/mm39/mm39"),
    split_name=False,
).with_columns(
    transcript_name=pl.col("transcript").apply(gtf_all.ts_to_gene),
)
# %%
t = res.filter(pl.col("name").str.starts_with("Abi3") & pl.col("aln_score").eq(60)).sort("aln_score")

# %%
res.filter(pl.col("name") == "Abi3_88").with_columns(
    transcript=pl.col("transcript").str.extract(r"(ENSMUST\d+)")
).join(fpkm, left_on="transcript", right_on="transcript_id(s)", how="left")
