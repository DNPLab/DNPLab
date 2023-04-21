"""Modules to process and analyse solid-state DNP eexperiments"""

import numpy as _np
from scipy import optimize
import dnplab as _dnp
import warnings
from scipy.constants import *

def process_DNP_profile(data, start, stop, frequency_sweep=True, normalize=False):

    """Cleanup DNPData object"""

    out = data.copy()

    if frequency_sweep == True:

        out.rename("t1", "freq")
        out.coords["freq"] = _np.linspace(start, stop, len(data.values))

    if frequency_sweep == False:

        out.rename("t1", "B0")
        out.coords["B0"] = _np.linspace(start, stop, len(data.values))

    if normalize == True:
        out = _dnp.normalize(out)

    out.attrs["experiment_type"] = "dnp_enhacement_profile"

    return out
