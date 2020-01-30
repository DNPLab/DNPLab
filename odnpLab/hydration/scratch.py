# hydration function(s)
# goal: take odnpImport and output tcorr and many observables

import numpy as np
from scipy import interpolate


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
    return interpolate.interp1d(power, T1,
                                kind='linear',
                                fill_value='extrapolate')


def getT1CorrectedEnhancement(E: np.array, power: np.array, t1p_fun: callable):
    """
    returns T1 corrected enhancements

    :param E: np.array of enhancements
    :param power: np.array of power in Watt unit, same length as E
    :param t1p_fun: callable function that returns T1 at any power
    :return: Ecorr, np.array of corrected enhancements with the same length as E
    """
    # TODO: Implement this
    raise NotImplementedError


def getKsigma(E: np.array, power: np.array):
    """
    returns k_sigma

    :param E: np.array of enhancements
    :param power: np.array of power in Watt unit, same length as E
    :return:
        float k_sigma
    """
    # TODO: Implement this
    raise NotImplementedError


def getTcorr(eps: float):
    """
    returns correlation time tcorr

    :param eps: float of epsilon
    :return:
        float tcorr
    """
    # TODO: Implement this
    raise NotImplementedError


def wip_getJH(ksi: float):
    """
    returns JH ksi_fit=((6*0.0063)-0.0063)/((6*0.0063)+(3*JH)+0.0063);
    :param ksi:
    :return:
    """
    # TODO: clean up
    raise NotImplementedError

