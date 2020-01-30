# hydration function(s)
# goal: take odnpImport and output tcorr and many observables

import numpy as np


def getT1p(T1: np.array, power: np.array):
    """
    returns a function to calculate T1 at arbitrary power

    :param T1: np.array of T1
    :param power: np.array of power in Watt unit, same length as T1
    :return
        a function that takes power and returns T1.

    """
    # TODO:// Implement this
    raise NotImplementedError


def getT1CorrectedEnhancement(E: np.array, power: np.array):
    """
    returns T1 corrected enhancements

    :param E: np.array of enhancements
    :param power: np.array of power in Watt unit, same length as E
    :return: Ecorr, np.array of corrected enhancements with the same length as E
    """
    # TODO:// Implement this
    raise NotImplementedError


def getEpsilon(E: np.array, power: np.array):
    """
    returns k_sigma

    :param E: np.array of enhancements
    :param power: np.array of power in Watt unit, same length as E
    :return:
        float k_sigma
    """
    # TODO:// Implement this
    raise NotImplementedError


def getTcorr(eps: float):
    """
    returns correlation time tcorr

    :param eps: float of epsilon
    :return:
        float tcorr
    """
    # TODO:// Implement this
    raise NotImplementedError


