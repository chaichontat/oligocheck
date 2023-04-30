import io

import pandas as pd

from oligocheck.merfish.alignment import run_bowtie
from oligocheck.merfish.external_data import ExternalData
from oligocheck.sequtils import parse_cigar, parse_sam, tm_hybrid


def calc_tm(row: pd.Series):
    seq = row["seq"][row["align_start"] : row["align_start"] + row["matches"]]
    return tm_hybrid(seq)


def check_align(ensembl: ExternalData, df: pd.DataFrame, gene: str | None = None):
    if gene is not None:
        tss = ensembl.get_transcripts(gene)
        check = df[df.gene == gene]
        check = check[~check.transcript.isin(tss)].copy()
    else:
        check = df.copy()

    check["align_start"], check["matches"] = zip(*check.cigar.map(parse_cigar))
    check["tm_offtarget"] = check.apply(calc_tm, axis=1)
    worst_case = (
        check.sort_values("matches").drop_duplicates("name", keep="last").sort_values(["tm_offtarget"])
    )
    # assert worst_case.iloc[40].tm < 40
    return worst_case


def get_blacklist(readout: str, ensembl: ExternalData):
    zeroth_readout = pd.read_csv(readout, index_col=0)
    fusedreadout = []
    for s1 in zeroth_readout.index:
        for s2 in zeroth_readout.index:
            fusedreadout.append(
                dict(
                    name=s1 + "_" + s2,
                    seq="AA" + zeroth_readout.loc[s1] + "AA" + zeroth_readout.loc[s2] + "AA",
                )
            )

    fusedreadout = pd.DataFrame(fusedreadout).set_index("name")["seq"]

    f = io.StringIO()
    for s1, row in zip(fusedreadout.index, fusedreadout):
        f.write(f"@{s1}\n")
        f.write(row.seq + "\n")
        f.write("+\n")
        f.write("~" * len(row.seq) + "\n")

    out = run_bowtie(f.getvalue(), "data/mm39/mm39", seed_length=16)
    f.close()

    zero = parse_sam(out, split_name=False).to_pandas()
    zero = zero[zero.transcript != "*"]
    zero = check_align(ensembl, zero)
    zero["bad_tm"] = zero.tm_offtarget > 33

    return set(
        map(
            tuple,
            (zero[zero.bad_tm].name.str.split("_", expand=True).to_numpy().tolist()),
        )
    )


# %%

# def check_align(df: pl.DataFrame, gene: str | None = None):
#     if gene is not None:
#         tss = ensembl.all_transcripts(gene)
#         check = df.filter(pl.col("gene") == "gene")
#         check = check.filter(pl.col("transcript").is_in(tss).is_not())
#     else:
#         check = df

#     check = check.with_columns(
#         [
#             pl.col("cigar")
#             .str.extract_all(r"(\d+)S(\d+)M")
#             .arr.eval(pl.element().cast(pl.UInt8))
#             .alias("matches")
#         ]
#     ).with_columns(pl.col("matches"))

#     # check["align_start"], check["matches"] = zip(*check.cigar.map(parse_cigar))
#     check["tm_offtarget"] = check.apply(calc_tm, axis=1)
#     worst_case = (
#         check.sort_values("matches").drop_duplicates("name", keep="last").sort_values(["tm_offtarget"])
#     )
#     # assert worst_case.iloc[40].tm < 40
#     return worst_case


# zero = parse_sam(out, split_name=False)
# zero = zero.filter(pl.col("transcript") != "*")
