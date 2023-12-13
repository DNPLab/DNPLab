# %%
"""
.. _plot_04_phase_cycling:

=======================================
Using the dnplab phase cycling function
=======================================

This example demonstrates how to use the phase cycling function on dnpdata objects.

"""
# The following example shows how and when the phase_cycle function can be conveniently used.
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

data = dnp.fourier_transform(data)["f2", (-900.0, -100.0)]

# %%
# we look at the real part of the spectra

data_real = data.real

# %%
# lets plot the spectrum for all 4 cycles

dnp.fancy_plot(data_real)

# %%
# clearly the spectra are phase cycled, but we'd like to have them phased all the same
# we can use dnp.phase_cycle for that:
# Note that [0(X),1(Y),2(-X),3(-Y)] gives the phase cycle Axes.

data_phased = dnp.phase_cycle(data, "Average", [0, 1, 2, 3])
dnp.fancy_plot(data_phased)
