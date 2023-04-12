# %%
import logging
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from functools import cache
from io import StringIO
from pathlib import Path
from typing import Iterable, Literal

import mygene
import pandas as pd
import primer3
import requests
from Bio import SeqIO

logger = logging.getLogger("nupack")
logger.propagate = False

from nupack import Model, SetSpec, Strand, Tube, tube_analysis
from nupack.analysis import Result, TubeResult

mg = mygene.MyGeneInfo()
seqs = SeqIO.index("/home/chaichontat/mer/mm39/Mus_musculus.GRCm39.cdna.all.fa", "fasta")


def reverse_complement(seq):
    return seq.translate(str.maketrans("ATCG", "TAGC"))[::-1]


@cache
def gen_model(
    t: float,
    formamide: float = 30,
    sodium: float = 0.3,
    magnesium: float = 0.0,
    **kwargs,
):
    return Model(
        material="dna",
        celsius=t + formamide * 0.65,
        sodium=sodium,
        magnesium=magnesium,
        **kwargs,
    )


def nonspecific_test(probe: str, seq: str, t: float = 37):
    model = gen_model(t)
    probe_ = Strand(reverse_complement("TTT" + probe + "TTT"), "probe")
    seq_ = Strand(seq, "seq")
    t1 = Tube(strands={seq_: 1e-10, probe_: 1e-9}, name="t1", complexes=SetSpec(max_size=2))
    # return tube_analysis(tubes=[t1], compute=["mfe"], model=model)
    result: Result = tube_analysis(tubes=[t1], compute=["mfe"], model=model)
    tube_result: TubeResult = result[t1]
    want = {}
    for complex, conc in tube_result.complex_concentrations.items():
        if complex.name == "(probe)":
            want["hairpin"] = conc / 1e-9
        elif complex.name == "(probe+seq)" or complex.name == "(seq+probe)":
            want["bound"] = conc / 1e-10
        elif complex.name == "(probe+probe)":
            want["dimer"] = conc / 1e-9
    return result, want


def hairpin(seq: str, t: float = 47):
    return nonspecific_test(seq, seq, t=t)


# %%
@cache
def get_seq(eid: str):
    if "." in eid:
        return str(seqs[eid].seq)  # type: ignore

    for i in range(20):
        try:
            return str(seqs[f"{eid}.{i}"].seq)  # type: ignore
        except KeyError:
            ...

    return requests.get(
        f"https://rest.ensembl.org/sequence/id/{eid}?type=cdna",
        headers={"Content-Type": "text/plain"},
    ).text


@cache
def gene_info(gene: str):
    link = f"https://rest.ensembl.org/lookup/symbol/mus_musculus/{gene}?expand=1"
    res = requests.get(link, headers={"Content-Type": "application/json"}).json()
    return res


@cache
def gene_to_eid(gene: str):
    return gene_info(gene)["id"]


@cache
def all_transcripts(gene: str):
    return [t["id"] for t in gene_info(gene)["Transcript"]]


@cache
def gene_to_transcript(gene: str):
    """Remove . from transcript first"""
    link = f"https://rest.ensembl.org/lookup/symbol/mus_musculus/{gene}?expand=1"
    res = requests.get(link, headers={"Content-Type": "application/json"}).json()
    return {
        "id": res["id"],
        "canonical": res["canonical_transcript"],
        "transcripts": [x["id"] for x in res["Transcript"]],
        "biotype": [x["biotype"] for x in res["Transcript"]],
    }


def transcripts_to_gene(ts: Iterable[str], species="mouse"):
    """Remove . from transcript first"""
    ts = [x.split(".")[0] for x in ts]
    res = mg.querymany(ts, fields="symbol", scopes="ensembl.transcript,symbol", species=species)
    return {x["query"]: x["symbol"] for x in res}


# %%
def tm(s: str) -> float:
    """Approximately the same as the results from the MATLAB script"""
    return primer3.calc_tm(s, mv_conc=300, dv_conc=0, dna_conc=1) + 5


formamide_molar = lambda percent: percent * 1000 * 1.13 / 45.04  # noqa: E731  # density / molar mass


@dataclass(frozen=True)
class Stringency:
    min_tm: float = 47
    max_tm: float = 52
    min_gc: float = 35
    max_gc: float = 65
    overlap: float = -1
    hairpin_tm: float = 30
    filter_quad_c: bool = True
    filter_quad_a: bool = True
    unique: bool = True  # to the extent that bowtie2 cannot detect


# fmt: off
@dataclass(frozen=True)
class Stringencies:
    # high   = Stringency(min_tm=49, min_gc=35, max_gc=65, hairpin_tm=40, filter_quad_c=True)
    high   = Stringency(min_tm=49, min_gc=35, max_gc=65, hairpin_tm=20, filter_quad_c=True,  unique=False)
    medium = Stringency(min_tm=49, min_gc=25, max_gc=70, hairpin_tm=30, filter_quad_c=True,  unique=False)
    low    = Stringency(min_tm=47, min_gc=35, max_gc=65, hairpin_tm=35, filter_quad_c=True, unique=False, overlap=15)
    low_ol = Stringency(min_tm=47, min_gc=25, max_gc=70, hairpin_tm=35, filter_quad_c=False, unique=False, overlap=15)
# fmt: on

# %%
if Path("gencode_vM32_transcripts.parquet").exists():
    everything = pd.read_parquet("gencode_vM32_transcripts.parquet")
else:
    gencode = pd.read_table(
        "/home/chaichontat/mer/gencode.vM32.primary_assembly.basic.annotation.gtf",
        comment="#",
        sep="\t",
        names=[
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ],
    )
    gencode = gencode[gencode.feature.isin(["transcript", "gene"])]
    gencode.attribute = gencode.attribute.apply(
        lambda row: dict([map(lambda t: t.replace('"', ""), kv.split(" ")) for kv in row.split("; ")])
    )
    trs_only = gencode[gencode.attribute.apply(lambda x: "transcript_id" in x)]
    everything = pd.DataFrame(trs_only.apply(lambda row: {**row, **row.attribute}, axis=1).to_list())
    everything = everything.drop(columns=["attribute"])
    everything.gene_id = everything.gene_id.map(lambda x: x.split(".")[0])
    everything.transcript_id = everything.transcript_id.map(lambda x: x.split(".")[0])
    everything.to_parquet("gencode_vM32_transcripts.parquet")


# gene = 'Pclaf'
# stg = 'high'
def block_bowtie(gene: str, tss: Iterable[str], stringency: Stringency, temp: Path):
    tss = [f"{gene}_{ts}" for ts in tss]

    for ts in tss:
        if not (temp / f"{ts}.fasta").exists():
            (temp / f"{ts}.fasta").write_text(f">{ts}\n{get_seq(ts.split('_')[1])}")

    [(temp / f"{ts}.fastq").unlink(missing_ok=True) for ts in tss]

    with ThreadPoolExecutor() as executor:
        list(
            executor.map(
                lambda ts: subprocess.run(
                    [
                        "python",
                        "blockParse.py",
                        "-f",
                        (temp / f"{ts}.fasta").as_posix(),
                        "-t",
                        str(stringency.min_tm),
                        "-T",
                        str(stringency.max_tm),
                        "-F",
                        "30",
                        "-O",
                    ],
                    # + (
                    #     ["-S", str(-stringency.overlap)]
                    #     if stringency.overlap <= 0
                    #     else []
                    # )
                    # + (["-O"] if stringency.overlap > 0 else []),
                    check=True,
                ),
                tss,
            )
        )

        # Need some time between python and bowtie2 for files to flush to disk
    while not all((temp / f"{ts}.fastq").exists() for ts in tss):
        print("Waiting for files to flush to disk...")
        time.sleep(0.2)

    with ThreadPoolExecutor() as executor:
        executor.map(
            lambda ts: subprocess.run(
                f"bowtie2 -x ../mm39/mm39 {(temp / f'{ts}.fastq').as_posix()} -t --no-hd -k 2 --local -D 20 -R 3 -N 1 -L 20 -i C,4 --score-min G,1,4 -S {(temp / f'{ts}.sam').as_posix()}",
                shell=True,
                check=True,
            ),
            tss,
        )


def combine_transcripts(gene: str, tss: Iterable[str], temp: Path):
    grand = []
    for ts in [f"{gene}_{ts}" for ts in tss]:
        file_read = [
            ",".join(line.strip().split("\t")[:10]) for line in (temp / f"{ts}.sam").read_text().split("\n")
        ]

        y = pd.read_csv(
            StringIO("\n".join(file_read)),
            sep=",",
            header=None,
            names=[
                "name",
                "flag",
                "transcript",
                "pos",
                "mapq",
                "cigar",
                "rnext",
                "pnext",
                "tlen",
                "seq",
            ],
        )
        y.transcript = y.transcript.apply(lambda x: x.split(".")[0])
        if "*" in y.transcript.value_counts().index:
            raise ValueError("Unmapped reads found")
        grand.append(y)

    grand = pd.concat(grand).drop_duplicates(subset=["transcript", "seq"], keep="first")
    return grand


def filter_specifity(tss_all: Iterable[str], grand: pd.DataFrame, stringency: Stringency, threshold=0.05):
    grand["ok"] = False
    grand["nonspecific_binding"] = 0.0
    tsss = set(tss_all)

    grand["ok"] = False
    grand.loc[grand.transcript.isin(tsss), "ok"] = True

    if not stringency.unique:
        for idx, row in grand[~grand.transcript.isin(tsss)].iterrows():
            seq = get_seq(row.transcript)
            ns = nonspecific_test(
                row.seq,
                seq[max(0, row.pos - 3) : min(row.pos + len(row.seq) + 3, len(seq))],
                37,
            )[1]["bound"]
            grand.loc[idx, ["nonspecific_binding", "ok"]] = (ns, ns < threshold)  # type: ignore

    res = []
    for _, rows in grand.groupby("seq"):
        if rows["ok"].all():
            rows["n_mapped"] = rows[rows.transcript.isin(tsss)].shape[0]
            res.append(rows)

    return pd.concat(res).drop(columns=["ok"])


def filter_candidates(res: pd.DataFrame, stringency: Stringency, n_cand=300):
    if stringency.overlap > 0 or len(res) < n_cand:
        # Not dropping anything here since we don't have enough.
        return res

    curr = res.n_mapped.max()
    picked = pd.DataFrame()
    # shuffle to get uniform coverage in case we get to n_cand early.
    while len(picked) < n_cand:
        shuffled = res[res.n_mapped == curr].sample(frac=1)
        if len(shuffled) + len(picked) > n_cand:
            shuffled = shuffled.iloc[: n_cand - len(picked)]
        picked = pd.concat([picked, shuffled]) if len(picked) else shuffled
    return picked


def calc_thermo(picked: pd.DataFrame):
    return picked.assign(
        hp=picked["seq"].map(
            lambda x: primer3.calc_hairpin_tm(x, mv_conc=300, dv_conc=0, dntp_conc=0, dna_conc=1)
            - 0.65 * 30  # formamide
        ),
        homodimer=picked["seq"].map(
            lambda x: primer3.calc_homodimer_tm(x, mv_conc=300, dv_conc=0, dntp_conc=0, dna_conc=1)
            - 0.65 * 30  # formamide
        ),
        tm=picked["seq"].map(
            lambda x: primer3.calc_tm(
                x,
                mv_conc=300,
                dv_conc=0,
                dna_conc=1,
                formamide_conc=formamide_molar(0.3),
            )
        ),
    )


def filter_thermo(picked: pd.DataFrame, stringency: Stringency):
    return picked[(picked.hp < stringency.hairpin_tm) & (picked.homodimer < 40)]


def handle_overlap(tss: Iterable[str], filtered: pd.DataFrame, stringency: Stringency):
    selected = set()
    for ts in tss:
        nool = sorted(
            [(r, idx, pos) for idx, r in filtered.iterrows() if (pos := get_seq(ts).find(r.seq)) > -1],
            key=lambda x: x[2],
        )
        curr = min(0, -stringency.overlap)
        if len(nool[1][0].seq) < stringency.overlap:
            raise ValueError("Overlap too large")
        for r, idx, pos in nool:
            if pos - curr >= 0:
                selected.add(idx)
                curr = pos + len(r.seq) - stringency.overlap
    return filtered.loc[list(selected)]


def cSpecStackFilter(seq):
    for i in range(6):
        if seq[i : i + 6].count("C") >= 4:
            return False
    return True


def ACompfilter(seq: str):
    return seq.count("T") / len(seq) < 0.28


# %%
# stg = "low"
# gene = "B2m"


def run(gene: str, stg: Literal["high", "medium", "low", "low_ol"]):
    stringency: Stringency = getattr(Stringencies, stg)

    eid = gene_to_eid(gene)
    tss_gencode = tuple(everything[everything.gene_id == eid]["transcript_id"].values)
    tss_all = all_transcripts(gene)

    print(f"Running {gene} with {len(tss_gencode)} transcripts")

    temp = Path("temp")
    temp.mkdir(exist_ok=True)

    block_bowtie(gene, tss_gencode, stringency, temp)
    grand = combine_transcripts(gene, tss_gencode, temp)
    res = filter_specifity(tss_all, grand, stringency)
    res = res.drop_duplicates(subset=["seq"], keep="first")
    # picked = filter_candidates(res, stringency)
    picked = calc_thermo(res)
    picked = filter_thermo(picked, stringency)
    # https://www.nature.com/articles/s41596-022-00750-2

    picked["ok_quad_c"] = ~picked["seq"].str.contains("GGGG")
    picked["ok_quad_a"] = ~picked["seq"].str.contains("TTTT")
    picked["ok_comp_a"] = picked["seq"].map(ACompfilter)
    picked["ok_stack_c"] = picked["seq"].map(cSpecStackFilter)

    m = handle_overlap(
        tss_gencode,
        picked,
        # picked[picked.ok_quad_c & picked.ok_quad_a & picked.ok_comp_a & picked.ok_stack_c],
        stringency,
    ).sort_values(by=["transcript", "pos"])
    return m


# %%


# %%


# return selected, filtered

# %%
# if __name__ == "__main__":
# selected, filtered = run("Pclaf", "high")


# %%


# picked = filter_candidates(res, stringency)
# picked = calc_thermo(picked, stringency)
# filtered = picked[(picked.hp < stringency.hairpin_tm) & (picked.homodimer < 55)]
# selected = handle_overlap(tss, filtered, stringency)

# if stringency.filter_quad_c:
#     selected = selected[~selected['seq'].str.contains("CCCC")]


# print(len(selected))
# %%

# -L ignore kmer < 2 c is counter bit. Want to be such that most kmers use only 1 counter, so 8.
# jellyfish count -m 18 -o output -t 32 -s 10G -L 2 -c 3 combi.fa
# jellyfish dump -c -L 2 output > output.txt
# count = dict()
# for i, x in jsorted.iterrows():
#     count[x[0]] = x[1]

# %%


# %%


# %%
# %%
# %%
# %%
# %%
# %%
# %%
