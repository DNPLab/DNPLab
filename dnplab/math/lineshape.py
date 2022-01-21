import numpy as np

def voigtian(x, x0, sigma, gamma):
    """
    Voigt is a combintaion of Gaussian and Lorentzian lineshapes
    """
    z = ((x0 - x) + 1j * gamma) / (sigma * np.sqrt(2.0))
    fit = np.real(wofz(z)) / (sigma * np.sqrt(2 * np.pi))
    return fit


def gaussian(x, x0, sigma):
    return (
        1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-((x - x0) ** 2) / (2 * sigma ** 2))
    )


def lorentzian(x, x0, gamma):
    return (1.0 / (np.pi * gamma)) * gamma ** 2 / ((x - x0) ** 2 + gamma ** 2)