# %%
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import click

from oligocheck.logging import log, setup_logging
from oligocheck.merfish.external_data import ExternalData, find_aliases

setup_logging()
# from prefect import flow, task
# from prefect.client import get_client


# async def set_limit():
#     async with get_client() as client:
#         await client.create_concurrency_limit(tag="parallel", concurrency_limit=8)
gtf = ExternalData(
    cache="data/mm39/gencode_vM32_transcripts.parquet",
    path="data/mm39/gencode.vM32.chr_patch_hapl_scaff.basic.annotation.gtf",
    fasta="data/mm39/combi.fa.gz",
)


# @task(tags=["parallel"], name="First pass", task_run_name="first-pass-{gene}")
def runpls(gene: str, skipdone: bool = False) -> str:
    # if not Path(f"output/{gene}_high.parquet").exists():
    if skipdone and Path(f"output/{gene}_final.parquet").exists():
        print(f"{gene} already done, skipping.")
        return ""

    log("running", gene)
    res = subprocess.run(
        ["python", "scripts/new_postprocess.py", gene],
        check=True,
        capture_output=True,
        cwd=os.getcwd(),
    )
    log(f"ran {gene}")
    return res.stdout.decode()


# %%


def check_gene_names(genes: list[str]):
    notfound = []
    ok: list[str] = []
    for gene in genes:
        try:
            gtf.gene_to_eid(gene)
            ok.append(gene)
        except ValueError:
            print(f"Gene {gene} not found in gtf")
            notfound.append(gene)
    converted, res = find_aliases(notfound)
    if len(res["dup"]) or len(res["missing"]):
        raise ValueError(f"Duplicated aliases {res['dup']} or missing aliases {res['missing']}")

    return ok + [x["symbol"] for x in converted.values()], {k: v["symbol"] for k, v in converted.items()}


# %%


@click.command()
@click.argument("genes_path", type=click.Path(exists=True))
@click.option("--skipdone", is_flag=True)
def main(genes_path: str | Path, skipdone: bool = False):
    # asyncio.run(set_limit())
    genes_path = Path(genes_path)
    genes: list[str] = list(filter(lambda x: x, Path(genes_path).read_text().splitlines()))
    sr = genes  # [gene for gene in genes if not Path(f"output/{gene}.parquet").exists()]
    sr, converted = check_gene_names(sr)
    if len(converted):
        (genes_path.parent / (genes_path.stem + "_convert.json")).write_text(json.dumps(converted))
        (genes_path.parent / (genes_path.stem + "_converted.txt")).write_text("\n".join(sorted(sr)))
    print(f"running {len(sr)} genes")
    print(sr)
    # run(genes)

    with ThreadPoolExecutor(32) as executor:
        for x in as_completed([executor.submit(lambda g: runpls(g, skipdone), gene) for gene in sr]):
            print("ok")
            try:
                x.result()
            except subprocess.CalledProcessError:
                print(x.result())


if __name__ == "__main__":
    main()
