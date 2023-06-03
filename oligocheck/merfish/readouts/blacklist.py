# %%

import polars as pl
import primer3

from oligocheck.merfish.alignment import gen_bowtie_index, gen_fasta, gen_fastq, run_bowtie
from oligocheck.sequtils import parse_sam

final = pl.read_csv("newreadouts.csv")
fusedreadout = []
for i in range(len(final)):
    for j in range(len(final)):
        fusedreadout.append(
            dict(
                name=f"{final[i, 'id']}_{final[j, 'id']}",
                seq="AA" + final[i, "seq"] + "AA" + final[j, "seq"] + "AA",
            )
        )

fused = pl.DataFrame(fusedreadout)
y = parse_sam(
    run_bowtie(
        gen_fastq(fused["name"], fused["seq"]).getvalue(),
        "data/humouse/humouse",
        seed_length=17,
        n_return=-1,
        threshold=19,
    ),
    split_name=False,
)
# %%
gen_bowtie_index(gen_fasta(fused["name"], fused["seq"]).getvalue(), "data/readouts", "newreadouts")

# %%
y = parse_sam(
    run_bowtie(
        gen_fastq(final["id"], final["seq"]).getvalue(),
        "data/readouts/newreadouts",
        seed_length=9,
        n_return=-1,
        threshold=15,
    ),
    split_name=False,
)
