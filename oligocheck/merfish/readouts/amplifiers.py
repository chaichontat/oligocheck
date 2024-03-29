# %%
from functools import cache
from itertools import cycle

import pandas as pd
import primer3
from Levenshtein import distance
from oligocheck.seqcalc import tm_q5

from oligocheck.sequtils import reverse_complement


@cache
def slide(x: str, n: int = 20) -> list[str]:
    return [x[i : i + n] for i in range(len(x) - n + 1)]


def min_dist(seqs: list[str], *, slider: str | None = None) -> int:
    if len(seqs) == 1:
        raise ValueError("Only one sequence provided")

    mini = min([distance(seqs[i], seqs[j]) for i in range(len(seqs)) for j in range(i + 1, len(seqs))])
    if slider is not None:
        for slid in slide(slider):
            for i, s in enumerate(seqs):
                if (d := distance(s, slid)) < mini:
                    print(slid, i, s, d)
                    mini = d

    return mini


# %% Branched amplification

primary_raw = """ACACTTTCACCTTCCCATTA TT ACACTTTCACCTTCCCATTA TT ACACTTTCACCTTCCCATTA TT ACACTTTCACCTTCCCATTA TT ACACTTTCACCTTCCCATTA TT ATCCTCCTTCAATACATCCC
TCCCAACACATCCTATCTCA TA TCCCAACACATCCTATCTCA TA TCCCAACACATCCTATCTCA TA TCCCAACACATCCTATCTCA TA TCCCAACACATCCTATCTCA TA ACACTACCACCATTTCCTAT
ACCCATTACTCCATTACCAT AT ACCCATTACTCCATTACCAT AT ACCCATTACTCCATTACCAT AT ACCCATTACTCCATTACCAT AT ACCCATTACTCCATTACCAT AT ACTCCACTACTACTCACTCT
TATCATCCTTACACCTCACT AA TATCATCCTTACACCTCACT AA TATCATCCTTACACCTCACT AA TATCATCCTTACACCTCACT AA TATCATCCTTACACCTCACT AA ACCCTCTAACTTCCATCACA
ACTCTCCTTTCCCATAACCT AT ACTCTCCTTTCCCATAACCT AT ACTCTCCTTTCCCATAACCT AT ACTCTCCTTTCCCATAACCT AT ACTCTCCTTTCCCATAACCT AT ACCACAACCCATTCCTTTCA
ATCATCAACACCTCATCAAC AT ATCATCAACACCTCATCAAC AT ATCATCAACACCTCATCAAC AT ATCATCAACACCTCATCAAC AT ATCATCAACACCTCATCAAC AT TTTCTACCACTAATCAACCC
ACCCTTTCAATCCATTCCAA TT ACCCTTTCAATCCATTCCAA TT ACCCTTTCAATCCATTCCAA TT ACCCTTTCAATCCATTCCAA TT ACCCTTTCAATCCATTCCAA TT ACCCTTTACAAACACACCCT
ACTTCTCACCTACCAATCAT AT ACTTCTCACCTACCAATCAT AT ACTTCTCACCTACCAATCAT AT ACTTCTCACCTACCAATCAT AT ACTTCTCACCTACCAATCAT AT TCCTATTCTCAACCTAACCT
TCCCACACATCATTTCCATT AA TCCCACACATCATTTCCATT AA TCCCACACATCATTTCCATT AA TCCCACACATCATTTCCATT AA TCCCACACATCATTTCCATT AA TATCCTTCAATCCCTCCACA
TCCTCCTATTCCCTAACAAC TA TCCTCCTATTCCCTAACAAC TA TCCTCCTATTCCCTAACAAC TA TCCTCCTATTCCCTAACAAC TA TCCTCCTATTCCCTAACAAC TA ACATTACACCTCATTCTCCC
TCCCATCTTCTTCCTAACTA TA TCCCATCTTCTTCCTAACTA TA TCCCATCTTCTTCCTAACTA TA TCCCATCTTCTTCCTAACTA TA TCCCATCTTCTTCCTAACTA TA TTTACTCCCTACACCTCCAA
TCACAATCCTACCACTACCT AA TCACAATCCTACCACTACCT AA TCACAATCCTACCACTACCT AA TCACAATCCTACCACTACCT AA TCACAATCCTACCACTACCT AA TTCTCCCTCTATCAACTCTA
TAACCTACACATCTCCACAA TA TAACCTACACATCTCCACAA TA TAACCTACACATCTCCACAA TA TAACCTACACATCTCCACAA TA TAACCTACACATCTCCACAA TA ACCCTTACTACTACATCATC
TCAATCCATTACCATCCCAC AT TCAATCCATTACCATCCCAC AT TCAATCCATTACCATCCCAC AT TCAATCCATTACCATCCCAC AT TCAATCCATTACCATCCCAC AT TCCTAACAACCAACTACTCC
AACCAAACCCACTACTACCA TT AACCAAACCCACTACTACCA TT AACCAAACCCACTACTACCA TT AACCAAACCCACTACTACCA TT AACCAAACCCACTACTACCA TT TCTATCATTACCCTCCTCCT
TACCCAACTAAACCCAACTC TA TACCCAACTAAACCCAACTC TA TACCCAACTAAACCCAACTC TA TACCCAACTAAACCCAACTC TA TACCCAACTAAACCCAACTC TA TATTCACCTTACAAACCCTC""".splitlines()

secondary_raw = """ATGAGGAAAGTGGTGTGAGA TT ATGAGGAAAGTGGTGTGAGA TT ATGAGGAAAGTGGTGTGAGA TT ATGAGGAAAGTGGTGTGAGA TT ATGAGGAAAGTGGTGTGAGA TA TAATGGGAAGGTGAAAGTGT
GAGGAGTGGATAAATGGTGT AT GAGGAGTGGATAAATGGTGT AT GAGGAGTGGATAAATGGTGT AT GAGGAGTGGATAAATGGTGT AT GAGGAGTGGATAAATGGTGT AT TGAGATAGGATGTGTTGGGA
AGTGTGGGATTGATGAGATA TT AGTGTGGGATTGATGAGATA TT AGTGTGGGATTGATGAGATA TT AGTGTGGGATTGATGAGATA TT AGTGTGGGATTGATGAGATA TT ATGGTAATGGAGTAATGGGT
TGTGGTTTGGAGATGATAGA TA TGTGGTTTGGAGATGATAGA TA TGTGGTTTGGAGATGATAGA TA TGTGGTTTGGAGATGATAGA TA TGTGGTTTGGAGATGATAGA TA AGTGAGGTGTAAGGATGATA
GAGATTAGAGATGAGTTGGA TA GAGATTAGAGATGAGTTGGA TA GAGATTAGAGATGAGTTGGA TA GAGATTAGAGATGAGTTGGA TA GAGATTAGAGATGAGTTGGA TA AGGTTATGGGAAAGGAGAGT
AGTTGAGGTGGGAGAGTATT AT AGTTGAGGTGGGAGAGTATT AT AGTTGAGGTGGGAGAGTATT AT AGTTGAGGTGGGAGAGTATT AT AGTTGAGGTGGGAGAGTATT AT GTTGATGAGGTGTTGATGAT
GGGTAGTGGGAATGATTTAT AT GGGTAGTGGGAATGATTTAT AT GGGTAGTGGGAATGATTTAT AT GGGTAGTGGGAATGATTTAT AT GGGTAGTGGGAATGATTTAT AT TTGGAATGGATTGAAAGGGT
AGGTAATGAGTTAGAGGTGA TT AGGTAATGAGTTAGAGGTGA TT AGGTAATGAGTTAGAGGTGA TT AGGTAATGAGTTAGAGGTGA TT AGGTAATGAGTTAGAGGTGA TT ATGATTGGTAGGTGAGAAGT
GGGATGTGATTTGTTAGGAA TT GGGATGTGATTTGTTAGGAA TT GGGATGTGATTTGTTAGGAA TT GGGATGTGATTTGTTAGGAA TT GGGATGTGATTTGTTAGGAA TT AATGGAAATGATGTGTGGGA
GATGAAGATTGAGGGAAGAA TT GATGAAGATTGAGGGAAGAA TT GATGAAGATTGAGGGAAGAA TT GATGAAGATTGAGGGAAGAA TT GATGAAGATTGAGGGAAGAA TT GTTGTTAGGGAATAGGAGGA
GGGATTATGGGTTTGTAGTA TT GGGATTATGGGTTTGTAGTA TT GGGATTATGGGTTTGTAGTA TT GGGATTATGGGTTTGTAGTA TT GGGATTATGGGTTTGTAGTA TA TAGTTAGGAAGAAGATGGGA
TAGAGGGAGTAAGATGAGGA TA TAGAGGGAGTAAGATGAGGA TA TAGAGGGAGTAAGATGAGGA TA TAGAGGGAGTAAGATGAGGA TA TAGAGGGAGTAAGATGAGGA TA AGGTAGTGGTAGGATTGTGA
GTGAAGTGGAAGGTGAGATT AT GTGAAGTGGAAGGTGAGATT AT GTGAAGTGGAAGGTGAGATT AT GTGAAGTGGAAGGTGAGATT AT GTGAAGTGGAAGGTGAGATT AA TTGTGGAGATGTGTAGGTTA
GAATGGAGGGTTAGAGGTAA TT GAATGGAGGGTTAGAGGTAA TT GAATGGAGGGTTAGAGGTAA TT GAATGGAGGGTTAGAGGTAA TT GAATGGAGGGTTAGAGGTAA TT GTGGGATGGTAATGGATTGA
TGGGATAGTATGTGGAAAGT AA TGGGATAGTATGTGGAAAGT AA TGGGATAGTATGTGGAAAGT AA TGGGATAGTATGTGGAAAGT AA TGGGATAGTATGTGGAAAGT AA TGGTAGTAGTGGGTTTGGTT
AGTTGGGTATGGAGAAAGGT AT AGTTGGGTATGGAGAAAGGT AT AGTTGGGTATGGAGAAAGGT AT AGTTGGGTATGGAGAAAGGT AT AGTTGGGTATGGAGAAAGGT AT GAGTTGGGTTTAGTTGGGTA""".splitlines()

# Binds to earlier stage on 3'
zeroth = [p[-20:] for p in primary_raw]
primary = [p[:20] for p in primary_raw]
secondary = [reverse_complement(s[:20]) for s in secondary_raw]
assert all([reverse_complement(p) == s[-20:] for p, s in zip(primary, secondary_raw)])


# %%
# Filter out anything with Levenshtein distance < 6
good_primary, good_secondary = [primary[0]], [secondary[0]]
for i in range(1, len(primary)):
    if min_dist([*good_primary, *good_secondary, primary[i], secondary[i]]) > 5:
        good_primary.append(primary[i])
        good_secondary.append(secondary[i])

assert min_dist([*good_primary, *good_secondary]) > 5
assert len(good_primary) == len(good_secondary)
# %%
# Reserved
# Full length from 900 onwards.
# First polyA readout
# 01-04: already chosen. Labeled as 00-03 in v1.

conditions = dict(
    mv_conc=300,
    dv_conc=0,
    dntp_conc=0,
    dna_conc=3,
)
bridges = pd.read_csv("ps_bridges.tsv", sep="\t")

dTmarker = "TTACACTCCATCCACTCAA"
selected = [
    dTmarker,
    "ACTCCAATCACCAAATAACA",
    "TCAATTTCACCACAAACCTC",
    "ACTAATCCTCTCAACACCAA",
    "TCTTCACCAAACCACTAAAC",
]
generated = []
dist = 6
for seq in bridges.iloc[:900].itertuples():
    for cutted in slide(seq.seq):
        if distance(cutted, selected[-1]) < dist:
            continue
        if not (49 <= primer3.calc_tm(cutted, **conditions) <= 53):
            continue

        # Distance check
        for s in [*selected, *good_primary, *good_secondary]:
            if distance(s, cutted) < dist:
                break
        else:
            if (
                primer3.calc_hairpin_tm(
                    gen := "AA" + reverse_complement(cutted) + "TAAAAAAAAAAAAAAAAAAAAAA",
                    **conditions,
                )
                <= 0
            ):
                selected.append(cutted)
                generated.append(gen)

# %%
# For generating PCR handles etc.
# Can only use 5 distances here since we're comparing 20-mers.
ok = []
for seq in bridges.iloc[900:].itertuples():
    for cutted in slide(seq.seq, 30):
        breakout = False
        # Distance check
        for cutted_s in slide(cutted, 20):
            for s in [*selected, *good_primary, *good_secondary]:
                if distance(s, cutted_s) < 5:
                    breakout = True
                    break
            if breakout:
                break
        else:
            ok.append(cutted)
# %%
sorted(list(zip(ok, map(tm_q5, ok))), key=lambda x: x[1])
# %%
sp6 = reverse_complement("ACGTGACTGCTCC" + SP6)
primary_rev = "CAAACTAACCTCCTTCTTCCTCCTTCCA"
secondary_rev = reverse_complement("TCACATCACACCTCTATCCATTATCAACCAC")

assert min_dist([*selected, *good_primary, *good_secondary], slider=primary_rev) > 4
assert (
    min_dist(
        [*selected, *good_primary, *good_secondary],
        slider=reverse_complement(secondary_rev),
    )
    > 4
)

for s in slide(primary_rev):
    for t in slide(reverse_complement(secondary_rev)):
        assert distance(s, t) > 4


# %%


# %%
existing_amp = len(good_primary)


def gen_amp(base: str, seq: str, reps: int = 5, rc: bool = False) -> str:
    if rc:
        seq = reverse_complement(seq)
        base = reverse_complement(base)
    fillers = cycle(["AA", "AT", "TA", "TT"])
    out = ""
    for f, _ in zip(fillers, range(reps)):
        out += f"{seq} {f} "
    return out + base


n = 48

zeroth_readouts = pd.DataFrame([{"name": f"Readout-{i:03d}-v2", "seq": selected[i]} for i in range(1, n + 1)])

primaries = [gen_amp(selected[i], p) for i, p in enumerate(good_primary, 1)]
primaries += [gen_amp(selected[i], selected[i + 300]) for i in range(existing_amp + 1, n + 1)]
primaries = [primary_rev + " " + p + " " + sp6 for p in primaries]


secondaries = [gen_amp(good_primary[i], s, rc=True) for i, s in enumerate(good_secondary)]
secondaries += [
    gen_amp(selected[i + 300], selected[i + 600], rc=True) for i in range(existing_amp + 1, n + 1)
]
secondaries = [secondary_rev + " " + s + " " + sp6 for s in secondaries]


assert min_dist([p.split(" ")[1] for p in primaries]) > 5
for p, s in zip(primaries, secondaries):
    assert p.split(" ")[1] == reverse_complement(s.split(" ")[-2])

secondaries_readout = good_secondary + [selected[i + 600] for i in range(existing_amp + 1, n + 1)]

assert all(p.split(" ")[-2] == selected[i] for i, p in enumerate(primaries, 1))
assert all(reverse_complement(p.split(" ")[1]) == s.split(" ")[-2] for p, s in zip(primaries, secondaries))
assert all(reverse_complement(s.split(" ")[1]) == fr for s, fr in zip(secondaries, secondaries_readout))

df_p = pd.DataFrame([{"Pool name": "PrimaryAmplifiers2", "Sequence": p} for p in primaries])
df_s = pd.DataFrame([{"Pool name": "SecondaryAmplifiers2", "Sequence": s} for s in secondaries])
out = pd.concat([df_p, df_s])

out.to_excel("opools.xlsx", index=False)

# Check that header is correct.
assert all(["C" in "".join(p.split(" ")[1:-1]) for p in primaries])
assert all(["G" not in "".join(p.split(" ")[1:-1]) for p in primaries])
assert all(["C" not in "".join(s.split(" ")[1:-1]) for s in secondaries])
assert all(["G" in "".join(s.split(" ")[1:-1]) for s in secondaries])

assert all([len(p) - reverse_complement(p).find(SP6) > 100 for p in primaries])
assert all([len(s) - reverse_complement(s).find(SP6) > 100 for s in secondaries])


# %%
# %% Prep for IDT
def genidt(seq: str) -> str:
    return "/5AmMC6/" + seq[:-2] + "*A*A"


idted = [genidt(x) for x in generated]
# %%

proteins = ["lectin", "mouse", "rabbit", "rat"]
for i, (name, sel, conj) in enumerate(zip(proteins, reversed(selected), reversed(idted)), 1):
    print(f"MERConjRev-{i:04d}-{name}\t{conj}\t\t100nm\tSTD")
    print(f"ReadoutRev-{i:04d}\t/5AmMC6/{sel}\t\t100nm\tSTD")
print()

# %%


# %%


def well_pos():
    for row in "ABCDEFGH":
        for col in range(1, 13):
            yield f"{row}{col:02d}"


frs = pd.DataFrame(
    [
        {"Well Position": well, "Name": f"FinalReadout-{i:02d}", "Sequence": "/5AmMC6/" + seq}
        for i, (well, seq) in enumerate(zip(well_pos(), secondaries_readout[:24]), 1)
    ]
).to_excel("final_readouts.xlsx", index=False)
#
for i, fr in enumerate(secondaries_readout, 1):
    print(f"FinalReadout-{i:02d}\t/5AmMC6/{fr}\t100nm\tSTD")
print()

# %%
# %%
