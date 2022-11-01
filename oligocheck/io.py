import hashlib
import urllib.request


def get_whitelist() -> list[str]:
    url = "https://assets.ctfassets.net/an68im79xiti/2jfcmwttryKsMq4AQqQcIS/539b681f94d6d5e1f3815d1817ce5448/CG000193_Barcode_Whitelist_forCustom_Feature_Barcoding_conjugates_RevA.txt"
    with urllib.request.urlopen(url) as response:
        seqs = response.read()

    if hashlib.md5(seqs).hexdigest() != "527656f6ac2b34232faf15e5016b1991":
        raise ValueError("Barcode file has changed.")

    return seqs.decode().splitlines()
