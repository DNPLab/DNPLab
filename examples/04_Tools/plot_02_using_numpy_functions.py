# %%
"""
.. _plot_02_using_numpy_functions:

=======================================
Using common numpy functions on DNPData
=======================================

This example demonstrates how to use common numpy funtions on DNPData.

"""
# %%
# Get Magnetic Resonance Properties
# =================================
# DNPLab stores a dictionary called ``gmrProperties`` with magnetic resonance properties of all nuclei of the periodic table. The dictionary is modeled after the Matlab function *gmr* written, and implemented by |GMRFunctionMatlab|. For more details see the detailed documentation of the dictionary. The dictionary stores the following parameters:
#
# * Spin Quantum Number
# * Gyromagnetic Ratio (Hz/T)
# * Nuclear Quadrupole Moment (fm^2, 100 barns)
# * Isotope Natural Abundance (%)
# * Relative Sensitivity with Respect to 1H at same B_{0}
# 
# (for nuclei with I > 1/2), and some more parameters. This dictionary can be used to provide nuclear properties in any calculation, it is also used by the DNPLab function ``mr_properties``.
# 
# To get started, first, setup the python environment:

import dnplab as dnp
import numpy as np
import matplotlib.pyplot as plt

# %%
# Let's load some example data

data=dnp.load("../../data/prospa/water_phase_cycled/data.2d")

# %% lets just plot it

data=np.mean(data,axis='Average')

dnp.fancy_plot(data)
dnp.plt.show()





