# %%
"""
.. _plot_03_analyze_odnp_data:

=================
Analyze ODNP Data
=================

This example demonstrates how to use the hydration module to analyze ODNP data.
"""
# %%
# Imports
# -------
import dnplab
import numpy as np
import matplotlib.pyplot as plt

# %%
# Setup your data dictionary
# --------------------------

data = {
    "E_array": np.array(),  # numpy array of signal enhancements
    "E_powers": np.array(),  # numpy array of microwave power levels used to collect enhancements
    "T1_array": np.array(),  # numpy array of T1 measurements in seconds
    "T1_powers": np.array(),  # numpy array of microwave power levels used to collect T1s
    "T10": 2.0,  # T1 measured with microwave power = 0
    "T100": 2.5,  # T1 measured for sample without unpaired spin and microwave power = 0
    "spin_C": 100e-6,  # spin concentration in M
    "magnetic_field": 0.35,  # magnetic field in T
    "smax_model": "tethered",  # choice of smax model or direct input of smax value
    "interpolate_method": "second_order",  # choice of interpolation method
}

# %%
# Standard constants
# ------------------
# In general the constants used in the calculations are kept the same as those found in literature. You may update these values if you wish but this is rarely necessary and will make direct comparisons with most existing literature invalid.

standard_constants = {
    "ksigma_bulk": 95.4,  # bulk ksigma value
    "krho_bulk": 353.4,  # bulk krho value
    "klow_bulk": 366,  # bulk klow value
    "tcorr_bulk": 54e-12,  # bulk tcorr value in seconds
    "D_H2O": 2.3e-9,  # bulk water diffusivity
    "D_SL": 4.1e-10,  # diffusivity of spin probe in bulk water
}

# %%
# Constants used in 2nd order T1 interpolation
# --------------------------------------------
# It is typically unnecessary to adjust these constants. Explanation can be found in the SI of https://doi.org/10.1021/jacs.1c11342

interpolation_constants = {
    "delta_T1_water": 1,  # change in water proton T1 due to microwaves
    "T1_water": 2.5,  # T1 of bulk water protons
    "macro_C": 100e-6,  # concentration of macromolecule in M
}

# %%
# combine the dictionaries of standard and additional constants into one,
constants = {**standard_constants, **interpolation_constants}

# %%
# Calculate results
# -----------------
# If any adjustments are made to the constants, pass both dictionaries to dnplab.hydration to return a dictionary of results,

results = dnplab.hydration(data, constants)

# %%
# If no adjustments are made to the constants you can skip the creation of the constants dictionary and pass just the data dictionary alone,

results = dnplab.hydration(data)

# %%
# Print a list of the calculated arrays and values,

print(results.keys())

# %%
# Plot the k_sigma array and fit, for example,

plt.figure()
plt.plot(data["E_powers"], results["ksigma_array"], "o", label="Data")
plt.plot(data["E_powers"], results["ksigma_fit"], label="Fit")
plt.legend()
plt.show()
