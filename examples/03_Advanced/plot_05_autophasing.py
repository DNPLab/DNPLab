# %%
"""
.. _plot_05_autophasing:

=======================================
Using the dnplab autophase function
=======================================

This example demonstrates how to use the autophase function on DNPData.

"""
# The following example shows how the autophase function can be conveniently used.
#
# To get started, first, setup the python environment:

import dnplab as dnp

# %%
# Let's load some example data
# the data consists of 4 fid that are phase cycled (0-90-180-270)

data = dnp.load("../../data/prospa/water_phase_cycled/data.2d")

# %%
# we are interested in the spectra
# and only the arbitrary part from -900 to -100 ppm

data = dnp.fourier_transform(data)['f2',(-900.0,-100.0)]

# %%
# we look at the real part of the spectra

data_real = data.real

# %%
# lets plot the spectrum for all 4 cycles

dnp.fancy_plot(data_real)

# %%
# clearly the spectra are phase cycled, but we'd like to have them phased all the same
# we can use dnp.autophase for that:

data_phased = dnp.autophase(data)
dnp.fancy_plot(data_phased)
