"""
01 - Import 1D Bruker Spectrum
==============================

This first, basic example demonstrates how to import a 1D FID into DNPLab and process the FID to obtain the NMR spectrum.

"""

# %%
# Script Header
# -------------
#
# Every DNPLab script will have at least one header row to import the DNPLab package. Importing DNPLab will also import other scripts (e.g. numpy, matplotlib, etc.) for internal use. You can also add other packages to be imported.

# sphinx_gallery_line_numbers = True
import sys
import numpy as np
sys.path.append("..")
import dnplab

# %%
# Import Data
# -----------
#
# Once the DNPLab package is imported, the FID can be loaded.

path = '../data/topspin/'
# folder = 20
data = dnplab.load(path, 'topspin', '20')


# %%
# Create Workspace
# ----------------
#
# A central concept of DNPLab is storing data in a workspace. In the following the workspace 'raw' and 'proc' is generated. The data object 'raw' contains the unprocessed spectral data. This is stored for reference, so the user can either choose to start working using the raw data or, once a few processing steps have been applied to the data, continue working on already processed data.

# Create workspace
ws = dnplab.create_workspace()
ws.add("raw", data)
ws.copy("raw", "proc")


# %%
# Process Data and Plot Spectrum
# ------------------------------
#
# First the FID is processed (remove DC offset, apply appodization and perform fourier transformation). Once the FID is process ti is plotted.

dnplab.dnpNMR.remove_offset(ws, {})
dnplab.dnpNMR.window(ws, linewidth = 10)
dnplab.dnpNMR.fourier_transform(ws, zero_fill_factor = 2)
dnplab.dnpNMR.autophase(ws, {})

dnplab.dnpResults.figure()
dnplab.dnpResults.plot(ws["proc"].real)
dnplab.dnpResults.xlim([50, -35])
dnplab.dnpResults.plt.xlabel("Chemical Shift [ppm]")
dnplab.dnpResults.plt.ylabel("Signal Amplitude [a.u.]")
dnplab.dnpResults.show()
