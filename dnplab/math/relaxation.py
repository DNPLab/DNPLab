import numpy as np

def t1_function(t, T1, M_0, M_inf):
    """Calculate exponential T1 curve

    .. math::
        f(t) = M_0 - M_{\infty} e^{-t/T_{1}}

    Args:
        t (array): time series
        T_{1} (float): T1 value
        M_{0} (float): see equation
        M_{\infty} (float): see equation

    Returns:
        array: T1 curve
    """

    return M_0 - M_inf * np.exp(-1.0 * t / T1)


def t2_function(t, M_0, T2, p):
    """Calculate stretched or un-stretched (p=1) exponential T2 curve

    .. math::
        f(t) = M_{0} e^{(-2(t/T_{2})^{p}}

    Args:
        t (array): time series
        M_{0} (float): see equation
        T_{2} (float): T2 value
        p (float): see equation

    Returns:
        array: T2 curve
    """

    return M_0 * np.exp(-2.0 * (t / T2) ** p)


def monoexp_fit(t, C1, C2, tau):
    """Calculate mono-exponential curve

    .. math::
        f(t) = C1 + C2 e^{-t/tau}

    Args:
        t (array): time series
        C1 (float): see equation
        C2 (float): see equation
        tau (float): see equation

    Returns:
        array: mono-exponential curve
    """

    return C1 + C2 * np.exp(-1.0 * t / tau)


def biexp_fit(t, C1, C2, tau1, C3, tau2):
    """Calculate bi-exponential curve

    .. math::
        f(t) = C1 + C2 e^{-t/tau1} + C3 e^{-t/tau2}

    Args:
        t (array): time series
        C1 (float): see equation
        C2 (float): see equation
        C3 (float): see equation
        tau1 (float): see equation
        tau2 (float): see equation

    Returns:
        array: bi-exponential curve
    """

    return C1 + C2 * np.exp(-1.0 * t / tau1) + C3 * np.exp(-1.0 * t / tau2)


def buildup_function(p, E_max, p_half):
    """Calculate asymptotic buildup curve

    .. math::
        f(p) = E_{max} * p / (p_{1/2} + p)

    Args:
        p (array): power series
        E_{max} (float): maximum enhancement
        p_{1/2} (float): power at half saturation

    Returns:
        array: buildup curve
    """

    return E_max * p / (p_half + p)