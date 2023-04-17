import fcntl
import json
from contextlib import contextmanager
from functools import cache
from pathlib import Path
from typing import Iterable

import mygene
import pandas as pd
import redis
import requests
from Bio import SeqIO

mg = mygene.MyGeneInfo()

seqs = SeqIO.index("/home/chaichontat/mer/oligocheck/data/mm39/Mus_musculus.GRCm39.cdna.all.fa", "fasta")

r_seq = redis.Redis(host="localhost", port=6379, db=0)
r_gene = redis.Redis(host="localhost", port=6379, db=1)


@contextmanager
def lock(path: str):
    locked_file_descriptor = open(Path(path), "w+")
    fcntl.lockf(locked_file_descriptor, fcntl.LOCK_EX)
    yield locked_file_descriptor
    locked_file_descriptor.close()


@cache
def get_seq(eid: str):
    try:
        eid, version = eid.split(".")
    except ValueError:
        version = None

    if (res := r_seq.get(eid)) is not None:
        return res.decode()

    if version is not None:
        try:
            r_seq.set(eid, res := str(seqs[f"{eid}.{version}"].seq))
            return res  # type: ignore
        except KeyError:
            ...

    for i in range(20):
        try:
            r_seq.set(eid, res := str(seqs[f"{eid}.{i}"].seq))  # type: ignore
            return res
        except KeyError:
            ...

    r_seq.set(
        eid,
        res := requests.get(
            f"https://rest.ensembl.org/sequence/id/{eid}?type=cdna",
            headers={"Content-Type": "text/plain"},
        ).text,
    )
    return res


gene_info_df = pd.DataFrame()


@cache
def gene_info(gene: str):
    if (res := r_gene.get(gene)) is not None:
        return json.loads(res.decode())

    link = f"https://rest.ensembl.org/lookup/symbol/mus_musculus/{gene}?expand=1"
    res = requests.get(link, headers={"Content-Type": "application/json"}).json()
    r_gene.set(gene, json.dumps(res))
    return res


@cache
def eid_info(eid: str):
    if (res := r_gene.get(eid)) is not None:
        test = json.loads(res.decode())
        if "error" not in test:
            return json.loads(res.decode())

    link = f"https://rest.ensembl.org/lookup/id/{eid}"
    res = requests.get(link, headers={"Content-Type": "application/json"}).json()
    if "error" not in res:
        r_gene.set(eid, json.dumps(res))
        return res
    raise ValueError(f"{res['error']} for {eid}")


@cache
def gene_to_eid(gene: str):
    return gene_info(gene)["id"]


@cache
def eid_to_gene(eid: str):
    if eid_info(eid)["biotype"] == "lncRNA":
        return "lncRNA"
    return eid_info(eid)["display_name"].split("-")[0]


@cache
def all_transcripts(gene: str):
    return [t["id"] for t in gene_info(gene)["Transcript"]]


@cache
def gene_to_transcript(gene: str):
    """Remove . from transcript first"""
    res = gene_info(gene)
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
            "/home/chaichontat/mer/oligocheck/data/mm39/gencode.vM32.chr_patch_hapl_scaff.basic.annotation.gtf",
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
    if not Path(path).exists():
        # Get tRNA and rRNA
        out = []
        ncrna = SeqIO.parse("data/mm39/Mus_musculus.GRCm39.ncrna.fa", "fasta")
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
    return set(pd.read_table(path, header=None)[0].str.split(">", expand=True)[1].dropna())
