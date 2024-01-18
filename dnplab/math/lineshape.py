import numpy as _np
from scipy.special import wofz
from ..constants import constants as _const

# import scipy.constants as _const


def voigtian(x, x0, sigma, gamma, integral=1.0, deriv=False):
    r"""Voigtian distribution. Lineshape given by a convolution of Gaussian and Lorentzian distributions.

    Args:
        x (array_like): input x
        x0 (float): center of distribution
        sigma (float): Gaussian Linewidth. Standard deviation of Gaussian distribution.
        gamma (float): Lorentzian linewidth. 2*gamma is the full width at half maximum (FWHM)
        integral (float): Integral of distribution
        deriv (boolean): Derivative of a Voigtian distribution (Gaussian broadened imaginary part of a phased spectrum).

    Returns:
        ndarray: Voigtian distribution

    The Voigtian distribution is defined as:

    .. math::

        f(x; x_0, \sigma, \gamma) = \frac{\operatorname{Re}[w(z)]}{\sigma \sqrt{2 \pi}}

    with

    .. math::
        z = \frac{x + i\gamma}{\sigma \sqrt{2}}

    Derivative:
    .. math::

        f(x) = \frac{1}{\sigma^3 \sqrt{2 \pi}} \left[ \gamma \operatorname{Im}[w(z)] - \left(x - x0\right) \operatorname{Re}[w(z)] \right]

    with

    .. math::
        z = \frac{\left( \left( x - x0 \right)  + 1j \gamma \right)}{\sigma \sqrt{2}}

    """
    if not isinstance(deriv, bool) and isinstance(deriv, (float, int)):
        deriv = True if int(deriv) == 1 else False
        # the scipy.optimize.curve_fit passing bool to an float

    if deriv == False:
        z = ((x0 - x) + 1j * gamma) / (sigma * _np.sqrt(2.0))
        out = _np.real(wofz(z)) / (sigma * _np.sqrt(2 * _const.pi))
        return integral * out

    elif deriv == True:
        z = ((x - x0) + 1j * gamma) / (sigma * _np.sqrt(2.0))
        xc = x - x0
        out = (
            1
            / sigma**3
            / _np.sqrt(2 * _const.pi)
            * (gamma * _np.imag(wofz(z)) - xc * _np.real(wofz(z)))
        )
        return integral * out

    else:
        raise ValueError("Derivative argument can only be either True or False.")


def gaussian(x, x0, sigma, integral=1.0):
    r"""Gaussian distribution.

    Args:
        x (array_like): input x
        x0 (float): Center of distribution
        sigma (float): Standard deviation of Gaussian distribution
        integral (float): Integral of distribution

    Returns:
        ndarray: Gaussian distribution

    The Gaussian distribution is defined as:

    .. math::

        f(x; x_0, \sigma) = \frac{1}{\sigma \sqrt{2 \pi}} \exp{\left(\frac{(x-x_0)^2}{2 \sigma^2}\right)}

    """
    return (
        integral
        / (sigma * _np.sqrt(2 * _const.pi))
        * _np.exp(-((x - x0) ** 2) / (2 * sigma**2))
    )


def lorentzian(x, x0, gamma, integral=1.0, deriv=False):
    r"""Lorentzian Distribution.

    Args:
        x (array_like): input x
        x0 (float): Center of distribution
        gamma (float): Lorentzian width. 2*gamma is full width at half maximum (FWHM)
        integral (float): Integral of distribution
        deriv (boolean): Derivative of a Lorentzian Distribution (Imaginary part of a phased spectrum)

    Returns:
        ndarray: Lorentzian distribution

    The Lorentzian distribution is defined as:

    .. math::

        f(x) = \frac{1}{\pi \gamma} \left[\frac{\gamma^2}{(x-x_0)^2 + \gamma^2}\right]

    Derivative:

    .. math::

        f(x) = \frac{1}{\pi \gamma} \left[\frac{- 2\gamma^2 (x-x_0)}{\left( (x-x_0)^2 + \gamma^2 \right)^2}\right]
    """
    if not isinstance(deriv, bool) and isinstance(deriv, (float, int)):
        deriv = True if int(deriv) == 1 else False
        # the scipy.optimize.curve_fit passing bool to an float

    if deriv == False:
        return (
            integral
            * (1.0 / (_const.pi * gamma))
            * gamma**2
            / ((x - x0) ** 2 + gamma**2)
        )

    elif deriv == True:
        return (
            integral
            * (-1.0 / (_const.pi * gamma))
            * gamma**2
            / ((x - x0) ** 2 + gamma**2) ** 2
            * 2.0
            * (x - x0)
        )

    else:
        raise ValueError("Derivative argument can only be either True or False.")
