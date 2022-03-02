import numpy as np
import warnings


class Lennard_Jones_Potential:
    def __init__(self) -> None:
        pass

    @staticmethod
    def potential(
        r: np.ndarray, 
        rc: float = None, 
        sigma: float = 1, 
        epsilon: float = 1,
        ignore_RuntimeWarning: bool = True
    ) -> np.ndarray:

        if ignore_RuntimeWarning:
            warnings.filterwarnings("ignore", category=RuntimeWarning)

        if type(r) is not np.ndarray:
            r = np.asarray(r, dtype='float64')

        s6 = (sigma/r)**6
        s12 = s6 * s6

        if rc is None:
            return 4*epsilon*(s12-s6)

        else:
            return np.where(
                np.logical_and(r<3, r!=0), 
                4*epsilon*((s12-s6) - ((sigma/rc)**12 - (sigma/rc)**6)), 
                0
            )
