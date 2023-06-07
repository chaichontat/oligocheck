# %%
from pathlib import Path

import polars as pl
from Levenshtein import distance

from oligocheck.merfish.alignment import gen_bowtie_index, gen_fasta, gen_fastq, run_bowtie
from oligocheck.merfish.filtration import count_match
from oligocheck.seqcalc import tm_fish
from oligocheck.sequtils import parse_sam

seed = "TTACACTCCATCCACTCAA"
df = pl.read_csv("data/readout_ref.csv", separator="\t")
# df = pl.read_csv("newreadouts.csv").with_columns(bit=pl.col("id"), name=pl.col("id").cast(pl.Utf8))

# %%
ok = [seed]
names = ["polyT"]
for row in df.iter_rows(named=True):
    bit, name, seq = row["bit"], row["name"], row["seq"]
    if distance(seq, seed) < 9:
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
offtarget = parse_sam(
    run_bowtie(
        gen_fastq(names=names, seqs=ok).getvalue(), "data/humouse/humouse", seed_length=12, threshold=14
    ),
    split_name=False,
)

# %%

offtarget = count_match(offtarget)

# %%
nottoobad = (
    offtarget.groupby("name").agg(pl.max("match_max")).sort("match_max").filter(pl.col("match_max") < 17)
)
picked2 = (
    picked.filter(pl.col("name").is_in(nottoobad["name"]))
    .with_columns(pl.col("seq").apply(tm_fish).alias("tm"))
    .sort("tm")
    .filter(pl.col("tm").is_between(51, 55))
)
picked2.sort("name").with_row_count("id", 1).write_csv("data/readout_ref_filtered.csv")
# %%

fusedreadout = []
for i in range(len(picked)):
    for j in range(len(picked)):
        fusedreadout.append(
            dict(
                name=f"{picked[i, 'name']}_{picked[j, 'name']}",
                seq="TT" + picked[i, "seq"] + "TT" + picked[j, "seq"] + "TT",
            )
        )

fused = pl.DataFrame(fusedreadout)
y = count_match(
    parse_sam(
        run_bowtie(
            gen_fastq(fused["name"], fused["seq"]).getvalue(),
            "data/humouse/humouse",
            seed_length=15,
            n_return=50,
            threshold=18,
        ),
        split_name=False,
    )
).with_columns(
    split1=pl.col("transcript").str.split("_").arr.get(0),
    split2=pl.col("transcript").str.split("_").arr.get(1),
)
# %%
y.groupby("name").agg(pl.col("match_max").max()).filter(pl.col("match_max") > 18).with_columns(
    split1=pl.col("name").str.split("_").arr.get(0),
    split2=pl.col("name").str.split("_").arr.get(1),
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
