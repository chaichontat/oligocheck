import hashlib
import shutil
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Literal, overload

import requests


@overload
def download(url: str, path: None = ..., return_bytes: Literal[True] = ...) -> bytes:
    ...


@overload
def download(url: str, path: None = ..., return_bytes: Literal[False] = ...) -> str:
    ...


@overload
def download(url: str, path: Path | str) -> None:
    ...


def download(url: str, path: str | Path | None = None, return_bytes: bool = False) -> str | bytes | None:
    with requests.get(url, stream=True) as r:
        if path is None:
            return r.content if return_bytes else r.text
        path = Path(path)
        with open(path / url.split("/")[-1], "wb") as f:
            shutil.copyfileobj(r.raw, f)


def get_whitelist() -> list[str]:
    url = "https://assets.ctfassets.net/an68im79xiti/2jfcmwttryKsMq4AQqQcIS/539b681f94d6d5e1f3815d1817ce5448/CG000193_Barcode_Whitelist_forCustom_Feature_Barcoding_conjugates_RevA.txt"
    seqs = download(url, return_bytes=True)
    if hashlib.md5(seqs).hexdigest() != "527656f6ac2b34232faf15e5016b1991":
        raise ValueError("Barcode file has changed.")

    return seqs.decode().splitlines()


def get_paintshop(path: Path | str = "data/paintshop", redownload: bool = False) -> None:
    urls = [
        "https://paintshop-bucket.s3.amazonaws.com/v1.2/resources/appending_plaintext/ps_if.txt",
        "https://paintshop-bucket.s3.amazonaws.com/v1.2/resources/appending_plaintext/ps_ir.txt",
        "https://paintshop-bucket.s3.amazonaws.com/v1.2/resources/appending_plaintext/ps_of.txt",
        "https://paintshop-bucket.s3.amazonaws.com/v1.2/resources/appending_plaintext/ps_or.txt",
    ]
    Path(path).mkdir(parents=True, exist_ok=True)
    if not redownload and len(list(Path(path).glob("*"))) == len(urls):
        return
    with ThreadPoolExecutor() as executor:
        executor.map(lambda url: download(url, path), urls)
