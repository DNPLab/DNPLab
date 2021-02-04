# %% [markdown]
"""
DNPLab basics tutorial
============================

This example demonstrates how to import and process data. The user also learns how to access the various components of the data object and generally how to interact with the functions.

"""
# %%
import numpy as np
import dnplab as dnp

# %%

# %% [markdown]
# First use some 1D topspin data. All data enters into the same object structure so
# this example applies to any NMR format.
data = dnp.dnpImport.load("../data/topspin/5")
# %%

# %% [markdown]
# Create a workspace, this is a place to store multiple objects if multiple processing
# methods or data manipulations are to be compared.
ws = dnp.create_workspace()
# %%

# %% [markdown]
# Add the data to the key 'raw'
ws.add("raw", data)
# %%

# %% [markdown]
# Copy the data to the processing buffer 'proc'
ws.copy("raw", "proc")
# %%

# %%
# Raw data are now in both ws["proc"] and ws["raw"]. The following processing steps will
# modify ws["proc"] but leave ws["raw"] untouched so that you can always return to the
# original data.
# %%

# %% [markdown]
# Perform some basic processing. Simply pass the entire workspace, 'ws', to the
# functions and the ws["proc"] will be manipulated. For example:
dnp.dnpNMR.remove_offset(ws)
dnp.dnpNMR.window(ws, linewidth=10)
dnp.dnpNMR.fourier_transform(ws, zero_fill_factor=2)
# %%

# %%
# This procedure removes the DC offset, multiplies by a window function, and fourier
# transforms.
# %%

# %% [markdown]
# You can manually phase correct using the 'autophase' function using
dnp.dnpNMR.autophase(ws, method="manual", phase=np.pi / 2)
# %%

# %% [markdown]
# Or use autophase to search for and apply the best phase angle
dnp.dnpNMR.autophase(ws)
# %%

# %% [markdown]
# If needed, access your processed spectrum as follows:
x_axis = ws["proc"].coords["f2"]  # ppm axis
spectrum = ws["proc"].values  # spectrum
# %%

# %% [markdown]
# By default, dimensions before FT are: direct is "t2", indirect is "t1".
# After FT, direct is "f2", indirect is still "t1". Rename dimensions at
# anytime with, for example,
ws["proc"].rename("f2", "ppm")
# %%

# %% [markdown]
# Access data attributes such as NMR frequency as follows:
nmr_frequency = ws["proc"].attrs["nmr_frequency"]
# %%

# %% [markdown]
# To see what attributes are available use
print(ws["proc"].attrs.keys())
# %%

# %% [markdown]
# Plot the frequency spectrum, imaginary and real
dnp.dnpResults.plot(ws["proc"].imag, label="Imaginary")
dnp.dnpResults.plot(ws["proc"].real, label="Real")
dnp.dnpResults.legend()
dnp.dnpResults.show()
# %%

# %% [markdown]
# Now suppose you have 2D inversion recovery data, import using
data = load("../data/topspin/304")
ws_2D = dnp.create_workspace()
ws_2D.add("raw", data)
ws_2D.copy("raw", "proc")
# %%

# %% [markdown]
# Process as above
dnp.dnpNMR.remove_offset(ws_2D)
dnp.dnpNMR.window(ws_2D, linewidth=10)
dnp.dnpNMR.fourier_transform(ws_2D, zero_fill_factor=2)
# %%

# %% [markdown]
# This time use a first order phase correction where the phase angle changes by linearly
# by delta=pi/16 over the spectral range and pivots 100 points into the spectrum.
dnp.dnpNMR.autophase(ws_2D, method="arctan", order="first", pivot=100, delta=np.pi / 16)
# %%

# %% [markdown]
# Now align the spectra
dnp.dnpNMR.align(ws_2D)
# %%

# %% [markdown]
# Store a copy of the data before they are reduced to integrals
ws_2D.copy("proc", "spec")
# %%

# %% [markdown]
# Now integrate
dnp.dnpTools.integrate(ws_2D)
# %%

# %% [markdown]
# By default the entire spectrum is integrated, be sure to optimize
# the center and width of the integration window. For example, first return
# to the data as it was before integrating by copying the saved 'spec'
# back to 'proc',
ws_2D.copy("spec", "proc")
# %%

# %% [markdown]
# Integrate over a 100ppm range centered at 10ppm
dnp.dnpTools.integrate(ws_2D, integrate_center=10, integrate_width=100)
# %%

# %% [markdown]
# Import the exponential fitting function with shortened syntax for convenience
from dnplab.dnpFit import exponential_fit as fit

# %%

# %% [markdown]
# Fit the integrals to an inversion recovery function to extract T1 using type="T1".
# Fit data are added to the object as the key 'fit'.
fit(ws_2D, type="T1")
# %%

# %% [markdown]
# If needed, access the fit data as follows:
x_axis = ws_2D["fit"].coords["t1"]  # indirect axis
fit_curve = ws_2D["fit"].values  # fit
T1 = ws_2D["fit"].attrs["T1"]  # T1 value
T1_stdd = ws_2D["fit"].attrs["T1_stdd"]  # standard deviation in T1
# %%

# %% [markdown]
# Plot the integrals and T1 fit
dnp.dnpResults.figure()
dnp.dnpResults.plot(ws_2D["proc"].imag, marker="o", linestyle="none", label="Imaginary")
dnp.dnpResults.plot(ws_2D["proc"].real, marker="o", linestyle="none", label="Real")
dnp.dnpResults.plot(ws_2D["fit"], label="Fit")
dnp.dnpResults.legend()
dnp.dnpResults.show()
# %%
