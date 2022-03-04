import numpy as np
import warnings
import numpy.typing as npt
from typing import Optional


class Lennard_Jones_Potential:
    def __init__(self) -> None:
        pass

    @staticmethod
    def potential(
        r_sqared: npt.ArrayLike, 
        rc: Optional[float] = None, 
        sigma: float = 1, 
        epsilon: float = 1,
        ignore_RuntimeWarning: bool = True
    ) -> np.ndarray:

        if ignore_RuntimeWarning:
            warnings.filterwarnings("ignore", category=RuntimeWarning)

        if type(r_sqared) is not np.ndarray:
            r_sqared = np.asarray(r_sqared, dtype='float64')

        s6 = (sigma*sigma/r_sqared)**3
        s12 = s6 * s6

        if rc is None:
            return 4*epsilon*(s12-s6)

        else:
            return np.where(
                np.logical_and(r_sqared<9, r_sqared!=0), 
                4*epsilon*((s12-s6) - ((sigma/rc)**12 - (sigma/rc)**6)), 
                0
            )
