"""Helper functions for numverical calculations"""

import numpy as _np
import dnplab as _dnp


def Jp(spin=1 / 2):
    """
    Single spin operator Jp (raising)

    Args:
        j:  Spin quantum number. Default is 1/2

    Returns:
        Jp: Jp single spin operator
    """

    multiplicity = int(2 * spin) + 1

    k = _np.linspace(-spin, spin - 1, multiplicity - 1)
    m = _np.sqrt(spin * (spin + 1) - k * (k + 1))
    out = _np.diag(m, 1)

    return out


def Jm(spin=1 / 2):
    """
    Single spin operator Jm (lowering)

    Args:
        j:  Spin quantum number. Default is 1/2

    Returns:
        Jm: Jm single spin operator
    """

    multiplicity = int(2 * spin) + 1

    k = _np.linspace(-spin, spin - 1, multiplicity - 1)
    m = _np.sqrt(spin * (spin + 1) - k * (k + 1))
    out = _np.diag(m, -1)

    return out


def Jz(spin=1 / 2):
    """
    Single spin operator Jz

    Args:
        j:  Spin quantum number. Default is 1/2

    Returns:
        Jz: Jz single spin operator
    """

    multiplicity = int(2 * spin) + 1
    m = _np.linspace(spin, -spin, multiplicity)
    out = _np.diag(m)

    return out


def Jx(spin=1 / 2):
    """
    Single spin operator Jx

    Args:
        j:  Spin quantum number. Default is 1/2

    Returns:
        Jx: Jx single spin operator
    """

    out = -0.5 * 1j * (Jp(spin) - Jm(spin))
    return out


def Jy(spin=1 / 2):
    """
    Single spin operator Jy

    Args:
        j:  Spin quantum number. Default is 1/2

    Returns:
        Jy: Jy single spin operator
    """

    out = 0.5 * (Jp(spin) + Jm(spin))
    return out
