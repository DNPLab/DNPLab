# %%
"""
.. _plot_05_indexing_dnpdata_objects:

==============================================
How to select a slice from a 2D dnpdata object
==============================================

This example demonstrates how to select a slice of a DNPData object.

You can index a DNPData object by specifying the name of the dimension and the index of the slice.
"""
# %%
# Import DNPLab and create a set of data
# --------------------------------------
# Use the lorentzian function to generate a 2d set of lorentzian distributions.

import numpy as np
from matplotlib.pylab import *

import dnplab as dnp

x = np.r_[-50:50:1024j]
y = np.r_[-10:10:1]

values = dnp.math.lineshape.lorentzian(x.reshape(-1, 1), y.reshape(1, -1), 0.5)

data = dnp.DNPData(values, ["f2", "sample"], [x, y])

# %%
# To specify a slice based on the index, we use an integer. This will select the slice at index 3.

data_slice_integer = data["sample", 3]  # get slice by index

# Taking the slice does not remove the sample dimension. We can remove dimensions of length 1 with the squeeze method
data_slice_integer.squeeze()  # remove "sample" dimension

# %%
# In many cases, we want the slice at a specific value of the coordinates. To do this, we use a float to specify the slice location. In python, by adding a period after the number, python interprets the number as a float instead of integer.

data_slice_float = data[
    "sample", 3.0
]  # get slice at index closest to coordinates value of 3.
data_slice_float.squeeze()  # again, we remove the "sample" dimension.

# %%
# Plot the result
# ---------------
# Let's plot the 1d slices:

figure()
dnp.plot(data_slice_integer)
dnp.plot(data_slice_float)
xlabel("Frequency (Hz)")
ylabel("Signal (a.u.)")
tight_layout()

# %%
# Similarly, we can also slice by specifying a range of values. To do this, we use a tuple specify the minimum and maximum values for the index.

data_slice_range = data["sample", (-3, 3)]

# %%
# For an advanced tutorial how indexing works and how to extract individual data slice from a multi-dimensional dnpData object see the :ref:`plot_02_extract_data` tutorial.
