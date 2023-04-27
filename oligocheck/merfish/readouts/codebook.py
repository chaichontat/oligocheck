# %%
import io
import subprocess
from pathlib import Path
from typing import Iterable, Sized

import numpy as np
import numpy.typing as npt


class CodebookPicker:
    def __init__(self, mhd4: npt.NDArray[np.bool_]) -> None:
        self.mhd4 = mhd4

    @staticmethod
    def _gen_codebook(mhd4: npt.NDArray[np.bool_], seed: int):
        rand = np.random.RandomState(seed)
        rmhd4 = mhd4.copy()
        rand.shuffle(rmhd4)
        return rmhd4

    def gen_codebook(self, seed: int):
        return self._gen_codebook(self.mhd4, seed)

    @staticmethod
    def _find_optimalish(mhd4: npt.NDArray[np.bool_], seed: int, fpkm: Sized):
        rmhd4 = CodebookPicker._gen_codebook(mhd4, seed)
        res = rmhd4[: len(fpkm)] * np.array(fpkm)[:, np.newaxis]
        tocalc = res.sum(axis=0)
        normed = tocalc / tocalc.sum()
        return -np.sum(normed * np.log2(normed)), tocalc

    def find_optimalish(self, fpkm: Sized, iterations: int = 5000):
        res = [self._find_optimalish(self.mhd4, i, fpkm)[0] for i in range(iterations)]
        best = np.argmax(res)
        return best, self._find_optimalish(self.mhd4, best, fpkm)[1]
