from functools import cache
from pathlib import Path
from typing import Iterable

import mygene
import pandas as pd
import requests
from Bio import SeqIO

mg = mygene.MyGeneInfo()
seqs = SeqIO.index("/home/chaichontat/mer/mm39/Mus_musculus.GRCm39.cdna.all.fa", "fasta")


@cache
def get_seq(eid: str):
    if "." in eid:
        try:
            return str(seqs[eid].seq)  # type: ignore
        except KeyError:
            eid = eid.split(".")[0]

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


def get_gencode(path: str = "gencode_vM32_transcripts.parquet"):
    if not Path(path).exists():
        gencode = pd.read_table(
            "/home/chaichontat/mer/mm39/gencode.vM32.chr_patch_hapl_scaff.basic.annotation.gtf",
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

    return pd.read_parquet(path)


def get_rrna(path: str) -> set[str]:
    if not Path("../mm39/rrna.fa").exists():
        # Get tRNA and rRNA
        out = []
        ncrna = SeqIO.parse("../mm39/Mus_musculus.GRCm39.ncrna.fa", "fasta")
        for line in ncrna:
            attrs = line.description.split(", ")[0].split(" ")
            actual = attrs[:7]
            description = " ".join(attrs[7:])
            attrs = [":".join(x.split(":")[1:]) if ":" in x else x for x in actual] + [
                str(line.seq),
                description,
            ]
            if len(attrs) < 8:
                attrs.append("")
            out.append(attrs)

        out = pd.DataFrame.from_records(
            out,
            columns=[
                "transcript_id",
                "type",
                "pos",
                "gene_id",
                "gene_biotype",
                "transcript_biotype",
                "gene_symbol",
                "seq",
                "description",
            ],
        )
        with open(path, "w") as f:
            for _, row in out[out.gene_biotype == "rRNA"].iterrows():
                f.write(f">{row.transcript_id}\n{row.seq}\n")

    return set(pd.read_table(path, header=None)[0].str.split(">", expand=True)[1].dropna())
