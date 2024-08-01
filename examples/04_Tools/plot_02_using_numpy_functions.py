# %%
"""
.. _plot_02_using_numpy_functions:

===============================================
Using common NumPy functions on dnpdata objects
===============================================

This example demonstrates how to use common numpy functions to manipulate dnpdata objects.

Note, this is still an experimental feature.

"""
# %%
# ===================================================
# How numpy array functions are operating on DNPData
# ===================================================
#
# Many numpy functions can directly be used on DNPData. What the function returns depends on the result:
#
# * when the result is a scalar a salar is returned
# * when the result is a ndarray a DNPData object is returned
#
# when the axis keyword is support by the numpy function one can provide the dimension (e.g. np.sum(mydata, axis='f2'))
# The corresponding axis is consumed and no longer in the returned DNPData object.
# The following example shows how this can be conveniently used.
#
# To get started, first, setup the python environment:

import dnplab as dnp
import numpy as np
import matplotlib.pyplot as plt

# %%
# Let's load some example data
# the data consists of 4 fid that are phase cycled (0-90-180-270)

data = dnp.load("../../data/prospa/water_phase_cycled/data.2d")

# %%
# we are interested in the spectra

data = dnp.fourier_transform(data)

# %%
# since we don't know what the spectra is made of, we want to have a quick look at the magnitude:

data_magn = np.abs(data)

# %%
# and lets plot the magnitude spectrum for all 4 cycles

dnp.fancy_plot(data_magn)

# %%
# since the spectra are phase cycled the mean of the real part of the spectrum should be 0.
# Let's check and plot that:

mean_real_spectrum = np.real(np.mean(data, axis="Average"))
total_mean = np.mean(mean_real_spectrum)
average_in_dims = "Average" in mean_real_spectrum.dims
print("The sum of the mean spectrum is {0} ".format(total_mean))
print("Average in mean_real_spectrum.dims: {0}".format(average_in_dims))
dnp.fancy_plot(mean_real_spectrum)

# %%
# the total mean could also be calculated directly using np.mean() without the axis keyword (up to numerical precision)
print(np.mean(np.real(data)), "~=", total_mean)
dnp.plt.show()
