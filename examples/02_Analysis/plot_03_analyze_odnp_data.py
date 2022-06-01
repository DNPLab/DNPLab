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
# Example Data
# ------------
enhancements = np.array(
    [
        0.57794113752189,
        -0.4688718613022250,
        -0.5464528159680670,
        -1.0725090541762200,
        -1.4141203961920700,
        -1.695789643686440,
        -1.771840068080760,
        -1.8420812985152700,
        -1.97571340381877,
        -2.091405209753480,
        -2.1860546327712800,
        -2.280712535872610,
        -2.4709892163826400,
        -2.5184316153191200,
        -2.556110148443770,
        -2.576413132701720,
        -2.675593912859120,
        -2.8153300703866400,
        -2.897475156648710,
        -3.0042154567120800,
        -3.087886507216510,
    ]
)

enhancement_powers = np.array(
    [
        0.0006454923080882520,
        0.004277023425898170,
        0.004719543572446050,
        0.00909714298712173,
        0.01344187403986090,
        0.01896059941058610,
        0.02101937603827090,
        0.022335737104727900,
        0.026029715703921800,
        0.02917012237740640,
        0.0338523245243911,
        0.03820738749745440,
        0.04733370907740660,
        0.05269608016472140,
        0.053790874615060400,
        0.05697639350179900,
        0.06435487925718170,
        0.07909179437004270,
        0.08958910066880800,
        0.1051813598911370,
        0.11617812912435900,
    ]
)

T1s = np.array(
    [
        2.020153734009,
        2.276836030132750,
        2.3708172489377400,
        2.4428968088189100,
        2.5709096032675700,
    ]
)

T1_powers = np.array(
    [
        0.000589495934876689,
        0.024242327290569100,
        0.054429505156431400,
        0.0862844940360515,
        0.11617812912435900,
    ]
)

# %%
# Setup your data dictionary
# --------------------------

data = {
    "E_array": enhancements,  # numpy array of signal enhancements
    "E_powers": enhancement_powers,  # numpy array of microwave power levels used to collect enhancements
    "T1_array": T1s,  # numpy array of T1 measurements in seconds
    "T1_powers": T1_powers,  # numpy array of microwave power levels used to collect T1s
    "T10": 2.0,  # T1 measured with microwave power = 0
    "T100": 2.5,  # T1 measured for sample without unpaired spin and microwave power = 0
    "spin_C": 100e-6,  # spin concentration in M
    "magnetic_field": 0.35,  # magnetic field in T
    "smax_model": "tethered",  # choice of smax model or direct input of smax value
    "interpolate_method": "second_order",  # choice of interpolation method
}

# %%
# Optional - adjust constants
# ---------------------------
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
# It is typically unnecessary to adjust the constants used in 2nd order interpolation but they are adjustable. Explanation can be found in the SI of https://doi.org/10.1021/jacs.1c11342

interpolation_constants = {
    "delta_T1_water": 1,  # change in water proton T1 due to microwaves
    "T1_water": 2.5,  # T1 of bulk water protons
    "macro_C": 100e-6,  # concentration of macromolecule in M
}

# %%
# Combine any dictionaries of constants into one with something like,

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
plt.grid()
plt.legend()
plt.show()
