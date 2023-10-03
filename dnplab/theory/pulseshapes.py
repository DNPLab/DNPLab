"""Functions to calculate AWG pulse shapes"""

import dnplab as _dnp
import numpy as _np
from scipy.linalg import expm as _expm

def chirp(tp, BW, B1, resolution=1.0e-9):
    """Calculate complex chirp pulse shape

    .. math::
        e^{i 2 \pi (k/2) (t - t_p/2)^2}

    Args:
        tp (float): Pulse length in (s)
        BW (float): Pulse excitation bandwidth in (Hz)

    Returns:
        tuple: tuple containing:
            t (*numpy.ndarray*): Time axis
            pulse (*numpy.ndarray*): Pulse shape
    """

    k = BW / tp
    t = _np.r_[0.0:tp:resolution]
    pulse = _np.exp(1.0j * 2.0 * _np.pi * ((k / 2.0) * ((t - tp / 2.0) ** 2.0)))

    pulse = pulse * B1

    return t, pulse


def wurst(tp, N, BW, B1, resolution=1.0e-9):
    """Calculate complex WURST pulse shape

    The WURST pulse can be constructed from two parts: 1) a hyperbolic secant, and 2) a frequency chirp.

    .. math::
        1 - \\text{abs} \left( \cos \left( \\frac{\pi}{t_p} (t - \\frac{t_p}{2}) + \\frac{\pi}{2} \\right) \\right) ^N

    Args:
        tp (float): Pulse length
        N (float): Exponent
        BW (float): Bandwidth in (Hz)
        B1 (float): Pulse strength in (Hz)
        resolution: Resolution (time increment) in (s). Default is 1 ns

    Returns:
        tuple: tuple containing:

        t (*numpy.ndarray*): Time axis
        pulse (*numpy.ndarray*): Real pulse shape of WURST pulse
    """

    t = _np.r_[0.0:tp:resolution]
    pulse_real = (
        1.0 - _np.abs(_np.cos(_np.pi * (t - tp / 2.0) / tp + _np.pi / 2.0)) ** N
    ) + 0j

    pulse_imag = (
        1.0 - _np.abs(_np.cos(_np.pi * (t - tp / 2.0) / tp + _np.pi / 2.0)) ** N
    ) + 0j

    # pulse_wurst = pulse_real + 1j * pulse_imag
    pulse_wurst = pulse_real

    # Use B1 = 1 here, otherwise the pulse shape is multiplied twice by B1
    t, pulse_chirp = _dnp.chirp(tp, BW, B1 = 1, resolution=resolution)

    pulse = B1 * pulse_wurst * pulse_chirp

    return t, pulse


def calculate_pulse_magnetization(t, pulse, BW, detection_op, spin = 1/2):
    """Evolve spin density matrix under pulse operator
    
    """

    M_list = []
    Mz_list = []

    npts = len(t)
    # dt = t[1] - t[0]
    dt = 1.e-9

    omega_array = _np.linspace(-BW/2, BW/2, npts)

    Mz_array = _np.zeros((len(t),len(omega_array)))

    for omega_ix,omega in enumerate(omega_array):
        
        sigma = _dnp.Jz(spin)

        Hz = 2*_np.pi * omega * _dnp.Jz(spin)
        
        for time_ix,time in enumerate(t):

            B1 = pulse[time_ix]
            H1 = _np.real(B1) * _dnp.Jx(spin) + _np.imag(B1) * _dnp.Jy(spin)

            H = Hz + H1

            P = _expm(1j*H*dt)                               # Define Propagator
            sigma = _np.dot(_np.dot(P,sigma),P.T.conj())      # Propagate Density Matrix

            Mz_value = _np.real(_np.trace(_np.dot(_dnp.Jz(spin),sigma)))

            Mz_array[time_ix, omega_ix] = Mz_value

        M = _np.trace(_np.dot(detection_op,sigma))                    # Detect

        Mz_value = _np.trace(_np.dot(_dnp.Jz(spin),sigma))

        M_list.append(M) # Append to FID array
        Mz_list.append(Mz_value)
    
    M = _np.array(M_list)
    Mz = _np.array(Mz_list)

    return M, Mz