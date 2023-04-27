# %%
import polars as pl
import primer3

from oligocheck.io import get_paintshop

pl.Config.set_fmt_str_lengths(100)
pl.Config.set_tbl_rows(100)
q5 = dict(dna_conc=500, dv_conc=5, dntp_conc=0.2, mv_conc=300)
paintshopprimers = ["ps_ir.txt", "ps_if.txt", "ps_or.txt", "ps_of.txt"]
get_paintshop()
pps = pl.concat([pl.read_csv(f"data/paintshop/{x}", separator="\t") for x in paintshopprimers])
# %%
tmss = []

candidates = pps.filter(pl.col("seq").str.contains(r"(CCC)|(GGG)").is_not())

candidates.with_columns(
    [
        pl.col("seq").apply(lambda x: primer3.calc_tm(x, **q5)).alias("tm"),
        pl.col("seq").apply(lambda x: primer3.calc_hairpin_tm(x, **q5)).alias("hairpin_tm"),
    ]
).sort("tm", descending=True).head(20)


# %% footer
(
    candidates.with_columns(
        [
            ("CCTATAGTGAGTCGTATTAG" + pl.col("seq").str.slice(0, 8)).alias("t7"),
            pl.col("seq").str.slice(0, 8).alias("t7chunk"),
        ]
    )
    .with_columns(
        [
            pl.col("t7").apply(lambda x: primer3.calc_tm(x, **q5)).alias("tm"),
            pl.col("t7").apply(lambda x: primer3.calc_hairpin_tm(x, **q5)).alias("hairpin_tm"),
        ]
    )
    .sort("tm", descending=True)
    .head(30)
)
# truncated.seq = truncated.seq.map(lambda x: x[:8])
# footers = []
# for i, (idx, p) in enumerate(truncated.iterrows()):
#     footers.append(
#         (
#             p["id"],
#             primer3.calc_tm("CCTATAGTGAGTCGTATTAG" + p.seq, **q5),
#             primer3.calc_hairpin_tm("CCTATAGTGAGTCGTATTAG" + p.seq, **q5),
#             "CCTATAGTGAGTCGTATTAG" + p.seq,
#         )
#     )
#     # for seq in out.seq:
#     #     tms.append(primer3.calc_hairpin_tm(reverse_complement(p) + seq[:30], **conditions))
#     # tmss.append((f"{i}r", max(tms)))

# footers = sorted(footers, key=lambda x: x[1])
# footers

# %%
# %%
