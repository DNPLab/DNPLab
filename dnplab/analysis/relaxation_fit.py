"""Modules which provide function to analyse relaxation measurements]"""

import numpy as _np
from ..fitting import fit
from ..math import *


def inversion_recovery_fit(integrals):
    # Estimate an initial guess from experimental data

    initial_guess = (2.0, -4000, 4000)

    fit_results = fit(relaxation.t1, integrals.real, dim="t1", p0=initial_guess)

    # fit returns dictionary with results

    print(fit_results)

    # print(fit_results['fit'])

    # fit = out['fit']
    # popt = out['popt']
    # err = out['err']

    # T1 = popt['popt',0]
    # M_0 = popt['popt',1]
    # M_inf = popt['popt',2]

    # return out
