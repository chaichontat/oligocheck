# %%
from concurrent.futures import ThreadPoolExecutor
from functools import reduce
from io import StringIO

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from tqdm import tqdm

from oligocheck.io import download


# %%
def get_publication(pub: str):
    return requests.get(
        f"https://www.encodeproject.org/publication-data/{pub}?frame=object",
        headers={"Accept": "application/json"},
    ).json()


x = get_publication("ENCSR574CRQ")


# %%
def get_encode(s: str):
    return requests.get(
        f"https://www.encodeproject.org/{s}?frame=object",
        headers={"Accept": "application/json"},
    ).json()


y = get_encode("/files/ENCFF919PTI/")

# %%
metadata = "https://www.encodeproject.org/documents/ab75e52f-64d9-4c39-aea0-15372479049d/@@download/attachment/ENCSR574CRQ_metadata.tsv"
metadata = requests.get(metadata).text
metadata = pd.read_csv(StringIO(metadata), sep="\t", index_col=0)

# %%

# %%
filtered = metadata[
    (metadata["File output type"] == "gene quantifications")
    & (metadata["Biosample term name"].str.find("brain") != -1)
]

with ThreadPoolExecutor(8) as executor:
    out = list(
        tqdm(
            executor.map(
                lambda x: get_encode(f"/files/{x}/"),
                filtered.index,
            ),
            total=len(filtered),
        )
    )
files_metadata = pd.DataFrame(out).set_index("accession")
files_metadata.to_csv("data/fpkm/metadata.tsv", sep="\t")

# %%
# download from azure uri
with ThreadPoolExecutor(8) as executor:
    out = list(
        tqdm(
            executor.map(
                lambda x: download(x, path="data/fpkm"),
                files_metadata["azure_uri"],
            ),
            total=len(files_metadata),
        ),
    )

# %%
files_metadata["simple_biosample_summary"].str.extract(r"(\d+\.?\d?) days")
# %%
useful_metadata = pd.DataFrame(
    dict(
        biosample=files_metadata["simple_biosample_summary"],
        age=files_metadata["simple_biosample_summary"].str.extract(r"(\d+\.?\d?) days")[0],
        tissue=filtered["Biosample term name"],
    ),
    index=files_metadata.index,
)
useful_metadata.to_csv("data/fpkm/useful_metadata.tsv", sep="\t")
# %%
dfs = [
    pd.DataFrame(pd.read_csv(f"data/fpkm/{x}.tsv", sep="\t", index_col="transcript_id(s)"))
    for x in useful_metadata[useful_metadata.age == "0"].index
]
assert all(df.index.is_unique for df in dfs)

# %%
combi = reduce(lambda x, y: x[["FPKM"]] + y[["FPKM"]], dfs)
combi.index = combi.index.str.split(".").str[0]  # remove version
combi.to_parquet("data/fpkm/P0_combi.parquet")
# %%
combi["FPKM"].hist(log=True, bins=np.logspace(0, 4, 100))
plt.gca().set_xscale("log")
# %%
# %%
