# %%
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import click
import pandas as pd

wants = [
    "Pclaf",
    "Cenpa",
    "Hells",
    "Top2a",
    "Mcm3",
    "Tipin",
    "Cenpf",
    "Mcm6",
    "Cdc20",
    "Ube2c",
    "Cdca8",
    "Birc5",
    "Smc2",
    "Rpa2",
    "Pcna",
    "Ccne1",
    "Tpx2",
    "Cks2",
    "Mki67",
    "Cdk1",
    "Prc1",
    "Lig1",
    "Ccnb1",
    "Mcm5",
    "Chaf1b",
    "Smc4",
    "Cenpe",
    "Dtl",
    "Uhrf1",
    "Spc25",
    "Nusap1",
    "Racgap1",
    "Chaf1a",
    "Mcm2",
    "Lmnb1",
    "Pbk",
    "Fen1",
    "Knstrn",
    "Mcm7",
    "Gmnn",
    "E2f1",
    "Tacc3",
    "Nasp",
    "Cdkn3",
    "Ccne2",
    "Clspn",
    "Ccnb2",
    "Incenp",
    "Ccng1",
    "Ezh2",
    "Topbp1",
    "Aurka",
    "Cenpq",
    "Kif22",
    "Aurkb",
    "Mcm4",
    "Mis18bp1",
    "Knl1",
    "Stmn1",
    "Ckap2",
    "Dbf4",
    "Psrc1",
    "Rif1",
    "Sgo2a",
    "Cdc6",
    "Fbxo5",
    "Rad21",
    "Ccnd3",
    "Kif4",
    "Sgo1",
    "Timeless",
    "Aspm",
    "Aunip",
    "Exo1",
    "Cenpw",
    "Ube2s",
    "Mad2l1",
    "Hjurp",
    "Kif11",
    "Nuf2",
    "Pmf1",
    "Cdkn2c",
    "Kif20b",
    "Esco2",
    "Cdc45",
    "Casp8ap2",
    "Nde1",
    "Cdca2",
    "Atad5",
    "Skp2",
    "Cdc25c",
    "Psmc3ip",
    "Tfdp1",
    "Cdkn1a",
    "Rpa1",
    "Ncapg",
    "Ckap5",
    "Ncapd3",
    "Bub3",
    "Kif2c",
    "Etaa1",
    "Vrk1",
    "Bub1b",
    "Cdk2",
    "Ncapg2",
]


def runpls(gene: str, stringency: str):
    # if not Path(f"output/{gene}_high.parquet").exists():
    # print("running", gene)
    subprocess.run(["genmer", gene, stringency, "--output", "output/"], check=True)


@click.command()
@click.argument("genes_path", type=click.Path(exists=True))
@click.argument("stringency")
def hello(genes_path: str, stringency: str):
    genes: list[str] = pd.read_csv(genes_path, header=None, dtype=str)[0]
    sr = [gene for gene in genes if not Path(f"output/{gene}_{stringency}.parquet").exists()]
    print(f"running {len(sr)} genes")

    with ThreadPoolExecutor(8) as executor:
        executor.map(runpls, sr, [stringency] * len(sr))


if __name__ == "__main__":
    hello()

# out = {}
# for gene in wants:
#     print(gene)
#     out[gene] = run(gene, "high")
# %%
# %%
# %%
