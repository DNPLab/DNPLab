import numpy as np
from scipy.special import wofz


def voigtian(x, x0, sigma, gamma, integral=1.0):
    r"""Voigtian distribution. Lineshape given by a convolution of Gaussian and Lorentzian distributions.

    :math:`f(x; x_0, \sigma, \gamma) = \frac{\operatorname{Re}[w(z)]}{\sigma \sqrt{2 \pi}}

    where,

    :math:`z = \frac{x + i\gamma}{\sigma \sqrt{2}}

    Args:
        x (array_like): input x
        x0 (float): center of distribution
        sigma (float): Gaussian Linewidth. Standard deviation of Gaussian distribution.
        gamma (float): Lorentzian linewidth. 2*gamma is the full width at half maximum (FWHM).
        integral (float): Integral of distribution

    Returns:
        ndarray: Voigtian distribution
    """
    z = ((x0 - x) + 1j * gamma) / (sigma * np.sqrt(2.0))
    out = np.real(wofz(z)) / (sigma * np.sqrt(2 * np.pi))
    return integral * out


def gaussian(x, x0, sigma, integral=1.0):
    r"""Gaussian distribution.

    :math:`f(x; x_0, \sigma) = \frac{1}{\sigma \sqrt{2 \pi}} \exp{\left(\frac{(x-x_0)^2}{2 \sigma^2}\right)}

    Args:
        x (array_like): input x
        x0 (float): Center of distribution
        sigma (float): Standard deviation of Gaussian distribution
        integral (float): Integral of distribution

    Returns:
        ndarray: Gaussian distribution
    """
    return (
        integral
        / (sigma * np.sqrt(2 * np.pi))
        * np.exp(-((x - x0) ** 2) / (2 * sigma ** 2))
    )


def lorentzian(x, x0, gamma, integral=1.0):
    r"""Lorentzian Distribution.

    f(x) = \frac{1}{\pi \gamma} \left[\frac{\gamma^2}{(x-x_0)^2 + \gamma^2}\right]

    Args:
        x (array_like): input x
        x0 (float): Center of distribution
        gamma (float): Lorentzian width. 2*gamma is full width at half maximum (FWHM).
        integral (float): Integral of distribution

    Returns:
        ndarray: Lorentzian distribution
    """
    return (
        integral * (1.0 / (np.pi * gamma)) * gamma ** 2 / ((x - x0) ** 2 + gamma ** 2)
    )
