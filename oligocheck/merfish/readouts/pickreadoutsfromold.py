# %%
from pathlib import Path

import polars as pl
from Levenshtein import distance

from oligocheck.geneframe import GeneFrame
from oligocheck.merfish.alignment import gen_bowtie_index, gen_fasta, gen_fastq, run_bowtie
from oligocheck.external.external_data import ExternalData
from oligocheck.seqcalc import tm_fish


seed = "TTACACTCCATCCACTCAA"
chosen = """RS0015
RS0083
RS0095
RS0332
RS0451
RS0468
RS0548
RS0639
RS0708
RS0763
RS0793
RS0820
RS0967
RS1040
RS1261
RS1470
RS3505
RS3811
RS4070
RS4111
RS4164
RS4216
RS4363
RS4465
RS0109
RS0584
RS0707
RS0730
RS1404
RS3791""".splitlines()
df = pl.read_csv("data/readout_ref.csv", separator="\t")
# df = pl.read_csv("newreadouts.csv").with_columns(bit=pl.col("id"), name=pl.col("id").cast(pl.Utf8))
gtf_all = ExternalData(
    cache="data/mm39/gencode_vM32_transcripts_all.parquet",
    path="data/mm39/Mus_musculus.GRCm39.109.gtf",
    fasta="data/mm39/combi.fa.gz",
)
fpkm = pl.read_parquet("data/fpkm/P0_combi.parquet")
# %%
ok = [seed, *df.filter(pl.col("name").is_in(chosen))["seq"]]
names = ["polyT", *df.filter(pl.col("name").is_in(chosen))["name"]]
for row in df.iter_rows(named=True):
    bit, name, seq = row["bit"], row["name"], row["seq"]
    if distance(seq, seed) < 7:
        continue
    for s in ok:
        if distance(s, seq) < 7:
            break
    else:
        ok.append(seq)
        names.append(name)
# %%
print(len(ok))
# %%
picked = pl.DataFrame({"name": names, "seq": ok})
# %%
offtarget = GeneFrame.from_sam(
    run_bowtie(
        gen_fastq(names=names, seqs=ok).getvalue(), "data/humouse/humouse", seed_length=12, threshold=14
    ),
    split_name=False,
)

# %%
nottoobad = (
    offtarget.groupby("name").agg(pl.max("match_max")).sort("match_max").filter(pl.col("match_max") < 17)
)
picked2 = (
    picked.filter(pl.col("name").is_in(nottoobad["name"]))
    .with_columns(pl.col("seq").apply(lambda x: tm_fish(x, formamide=0)).alias("tm"))
    .filter(pl.col("tm").is_between(51, 56.5))
    .sort("name")
)


# %%


ddf = pl.DataFrame({"name": chosen}).with_row_count("id", 1).join(df[["name", "seq"]], on="name").sort("id")
ddf.write_csv("data/readout_ref_filtered.csv")

ddf.with_columns(seq="/5AmMC6/" + pl.col("seq") + "/3AmMO/")[["name", "seq"]].write_csv(
    "readouts_idt.tsv", separator="\t"
)


# %%

fusedreadout = []
for i in range(len(picked)):
    # for j in range(len(picked)):
    fusedreadout.append(
        dict(
            name=f"{picked[i, 'name']}_{picked[i, 'name']}",
            seq="AA" + picked[i, "seq"] + "AA" + picked[i, "seq"] + "AA",
        )
    )

fused = pl.DataFrame(fusedreadout)
y = GeneFrame.from_sam(
    run_bowtie(
        gen_fastq(fused["name"], fused["seq"]).getvalue(),
        "data/mm39/mm39",
        seed_length=11,
        n_return=50,
        threshold=14,
    ),
    split_name=False,
).with_columns(
    split1=pl.col("transcript").str.split("_").list.get(0),
    split2=pl.col("transcript").str.split("_").list.get(1),
)
# %%

y.sort("match_max").filter(pl.col("match_max") > 16).filter(pl.col("flag") & 16 == 0).with_columns(
    transcript_name=pl.col("transcript").apply(gtf_all.ts_to_tsname)
)
# %%
y.groupby("name").agg(pl.col("match_max").max()).filter(pl.col("match_max") > 18).with_columns(
    split1=pl.col("name").str.split("_").list.get(0),
    split2=pl.col("name").str.split("_").list.get(1),
).write_csv("data/readout_fused_bad.csv")


# %%
# %%
# y.groupby("name").agg(pl.col("match_max").max()).to_pandas().hist()
gen_bowtie_index(gen_fasta(fused["name"], fused["seq"]).getvalue(), "data/readouts", "newreadouts")

# %%
y = count_match(
    parse_sam(
        run_bowtie(
            gen_fastq(picked["name"], picked["seq"]).getvalue(),
            "data/readouts/newreadouts",
            seed_length=9,
            n_return=-1,
            threshold=14,
        ),
        split_name=False,
    )
)

# %%


y.filter(pl.col("name").ne(pl.col("split1")) & pl.col("name").ne(pl.col("split2")))


# %%
def wells(n: int):
    assert n > 0
    col = (n - 1) // 8
    row = (n - 1) % 8
    return f"{chr(ord('A') + row)}{col + 1}"


df = (
    pl.read_csv("data/readout_ref_filtered.csv")
    .filter(pl.col("name").str.starts_with("RS"))
    .with_columns(pl.col("id").apply(wells).alias("Well Position"))
    .with_columns(seq="/5AmMC6/" + pl.col("seq") + "/3AmMO/")
).rename({"name": "Name", "seq": "Sequence"})[["Well Position", "Name", "Sequence"]]
df.write_excel("data/readout_ref_filtered.xlsx")
# %%
