# %%
import pandas as pd
import primer3

from oligocheck.io import get_paintshop

q5 = dict(dna_conc=500, dv_conc=5, dntp_conc=0.2, mv_conc=300)
# tosearch = combi.seq.map(lambda x: x[:30]).unique()
paintshopprimers = ["ps_ir.txt", "ps_if.txt", "ps_or.txt", "ps_of.txt"]
get_paintshop()
pps = pd.concat([pd.read_csv(f"data/paintshop/{x}", sep="\t") for x in paintshopprimers])
tmss = []
for i, (idx, p) in enumerate(pps.iterrows()):
    if "CCC" in p.seq or "GGG" in p.seq:
        continue
    tmss.append((p["id"], primer3.calc_tm(p.seq, **q5), primer3.calc_hairpin_tm(p.seq, **q5), p.seq))
    # for seq in out.seq:
    #     tms.append(primer3.calc_hairpin_tm(reverse_complement(p) + seq[:30], **conditions))
    # tmss.append((f"{i}r", max(tms)))

tmss = sorted(tmss, key=lambda x: x[1])
tmss
# %% footer
truncated = pps.copy()
truncated.seq = truncated.seq.map(lambda x: x[:8])
footers = []
for i, (idx, p) in enumerate(truncated.iterrows()):
    footers.append(
        (
            p["id"],
            primer3.calc_tm("CCTATAGTGAGTCGTATTAG" + p.seq, **q5),
            primer3.calc_hairpin_tm("CCTATAGTGAGTCGTATTAG" + p.seq, **q5),
            "CCTATAGTGAGTCGTATTAG" + p.seq,
        )
    )
    # for seq in out.seq:
    #     tms.append(primer3.calc_hairpin_tm(reverse_complement(p) + seq[:30], **conditions))
    # tmss.append((f"{i}r", max(tms)))

footers = sorted(footers, key=lambda x: x[1])
footers

# %%
