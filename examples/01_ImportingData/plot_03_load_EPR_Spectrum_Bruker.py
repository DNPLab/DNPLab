# %%
"""
.. _plot_03_load_EPR_Spectrum_Bruker:

======================================
Load EPR spectrum in Bruker EMX format
======================================

In this example we demonstrate how to load and EPR spectrum and process the data.

"""
# %%
# Load EPR Data
# -------------
# First, add DNPLab to the Python environment,

import dnplab as dnp

# %%
# and then import an EPR spectrum. DNPLab can handle spectra recorded on different spectrometers such as the Bruker ElexSys, the Bruker EMX system, or home-built spectrometers running on Boris Epel's software SpecMan4EPR. In this example we will load a spectrum from a Bruker EMX system.

data = dnp.load("../../data/bes3t/1D_CW.DTA")

# %%
# Process EPR Data
# ----------------
# In this section, we will demonstrate some basic EPR processing.

# %%
# First, let's perform a baseline correction using a zeroth order polynomial to remove a DC offset:
data_proc = dnp.remove_background(data, dim="B0")

# %%
# Here a new dnpData object is created containing the corrected data. This is helpful, if the processing for different data sets need to be compared. The remove_background function will calculate a zero order polynomial background and will subtract this value from the data. To plot the corrected spectrum simply use:

dnp.fancy_plot(data_proc, xlim=[344, 354], title="EPR Spectrum")

# %%
# The ''fancy_plot'' function is very helpful to create simple plots. For more complicated figures the matplotlib functions can be used. Note, that the plotting functions of the matplotlib package are already loaded into the DNPLab environment.

dnp.plt.figure()
dnp.plt.plot(data.coords["B0"], data.values.real, label="No Background Correction")
dnp.plt.plot(
    data_proc.coords["B0"], data_proc.values.real, label="Background Correction"
)
dnp.plt.xlabel("Magnetic Field (mT)")
dnp.plt.ylabel("EPR Signal Intensity (a.u.)")
dnp.plt.grid(True)
dnp.plt.tight_layout()
dnp.plt.legend()
dnp.plt.show()

# %%
# Note the DC offset of about -0.5.

# %%
# Show EPR Attributes
# -------------------
# To show a list of attributes with the EPR spectrum

dnp.fancy_plot(data_proc, xlim=[344, 354], title="EPR Spectrum", showPar=True)
dnp.plt.show()
