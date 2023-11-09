from scipy.constants import mu_0, c, pi, epsilon_0, hbar, h, N_A, m_e, eV
from scipy.constants import e as e_charge

from scipy.constants import physical_constants as pc

e_gyro = pc["electron gyromag. ratio in MHz/T"][0]
p_gyro = pc["proton gyromag. ratio in MHz/T"][0]
mub = pc["Bohr magneton"][0]
mub_Hz = pc["Bohr magneton in Hz/T"][0]


__all__ = [
    "mu_0",
    "c",
    "pi",
    "epsilon_0",
    "hbar",
    "h",
    "N_A",
    "m_e",
    "e_charge",
    "e_gyro",
    "p_gyro",
    "eV",
    "mub",
    "mub_Hz",
]
