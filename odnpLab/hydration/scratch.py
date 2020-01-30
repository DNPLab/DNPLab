# hydration function(s)
# goal: take odnpImport and output tcorr and many observables

import numpy as np
from scipy import interpolate
from scipy import optimize


def getT1p(T1: np.array, power: np.array):
    """
    returns a function to calculate T1 at arbitrary power

    :param T1: np.array of T1
    :param power: np.array of power in Watt unit, same length as T1
    :return
        a function that takes power and returns T1.
        points inside the data range will be linearly interpolated
        points outside the data range will be extrapolated
    """
    assert len(T1) == len(power)
    return interpolate.interp1d(power, T1,
                                kind='linear',
                                fill_value='extrapolate')


def getKsigma(E: np.array, power: np.array, t1p_fun: callable):
    """
    returns k_sigma

    :param E: np.array of enhancements
    :param power: np.array of power in Watt unit, same length as E
    :param t1p_fun: callable function that returns T1 at any power
    :return:
        float k_sigma
    """
    # TODO: Implement this
    raise NotImplementedError


def getTcorr(ksi: float):
    """
    returns correlation time tcorr

    :param ksi: float of epsilon
    :return:
        float tcorr
    """
    def get_ksi(tcorr:float):
        """
        returns ksi for any given tcorr
        :param tcorr: float
        :return:
            float ksi
        """
        # TODO: Implement this
        raise NotImplementedError

    # root finding
    # see https://docs.scipy.org/doc/scipy/reference/optimize.html
    results = optimize.root_scalar(
        lambda tcorr: sum((get_ksi(tcorr) - ksi) ** 2),
        method='newton',
        x0=0)
    assert results.converged
    return results.root
