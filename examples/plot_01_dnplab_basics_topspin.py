# %%
"""
.. _01_dnplab_baseics_topspin_data:

=================================
01 - DNPLab Basics (TopSpin Data)
=================================

This example demonstrates how to import and process data recorded using TopSpin. The user also learns how to access the various components of the data object and generally how to interact with the functions.

"""
# %%
# To use DNPLab you need to import dnplab. In this example, also numpy is imported. However, this is only necessary if you would like to access numpy functions and procedures directly.

import dnplab as dnp
import numpy as np

# %%
# Load and Process the NMR Spectrum
# ---------------------------------
# First use some sample 1D TopSpin data. Example data is located in the data folder. All data enters into the same object structure so this example applies to any NMR format.

data = dnp.dnpImport.load("../data/topspin/5")

# %%
# Next, create a workspace, this is a place to store multiple objects (e.g. copies of your data set) if multiple processing methods or data manipulations are to be compared.

ws = dnp.create_workspace()

# %%
# This will create a workspace variable called ws, and is handled as a Python dictionary. Next, add the data to the key 'raw'

ws.add("raw", data)

# %%
# Copy the data to the processing buffer 'proc'

ws.copy("raw", "proc")

# %%
# Raw data are now in both ws["proc"] and ws["raw"]. The following processing steps will only modify ws["proc"] but leaves ws["raw"] untouched. That way, you can always return to the original data. Furthermore, when you save the workspace (ws), you will always safe a copy of the raw data.
#
# Next, perform some basic processing. Simply pass the entire workspace, 'ws', to the functions and ws["proc"] will be manipulated. For example:

dnp.dnpNMR.remove_offset(ws)
dnp.dnpNMR.window(ws, linewidth = 20)
dnp.dnpNMR.fourier_transform(ws, zero_fill_factor = 2)

# %%
# This procedure removes the DC offset, performs a apodization with a window function (20 Hz linewidth), followed by a Fourier transforms.
#
# NMR spectra are typically phase sensitive and the final processing step is often phasing the spectra. In DNPLab you can either manually phase correct the spectrum using the 'autophase' function using the 'manual' keyword and giving a phase explicitely:

dnp.dnpNMR.autophase(ws, method = "manual", phase = np.pi / 2)

# %%
# Or by using the autophase function to automatically search for and apply the best phase angle

dnp.dnpNMR.autophase(ws)

# %%
# Plot Processed Data
# -------------------
# If needed, access your processed spectrum as follows:

x_axis = ws["proc"].coords["f2"]        # ppm axis
spectrum = ws["proc"].values            # spectrum

# %%
# DNPLab automatically takes care of the name of the coordinates. By default, the dimensions are named "t2" for the direct dimensions and "t1" for the indirect dimension before performing a Fourier Transformation. After the Fourier Transformation these dimensions will be called "f2" and "f1". However, you can rename dimensions at anytime. For example, like this:

ws["proc"].rename("f2", "ppm")

# %%
# Access data attributes such as the NMR frequency as follows:

nmr_frequency = ws["proc"].attrs["nmr_frequency"]

# %%
# The following will generate a plot of the NMR spectrum, and will show the imaginary and real part in the same figure.

dnp.dnpResults.figure()
dnp.dnpResults.plot(ws["proc"].imag, label="Imaginary")
dnp.dnpResults.plot(ws["proc"].real, label="Real")
dnp.dnpResults.xlim([100, -100])
dnp.dnpResults.plt.xlabel("Chemical Shift [ppm]")
dnp.dnpResults.plt.ylabel("Signal Amplitude [a.u.]")
dnp.dnpResults.legend()
dnp.dnpResults.plt.grid(True)
dnp.dnpResults.show()
