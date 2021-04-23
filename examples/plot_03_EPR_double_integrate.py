# %%
"""
.. _03_dnplab_epr_tutorial:

========================
03 - DNPLab EPR Tutorial
========================

This example demonstrates how to use DNPLab to load an EPR spectrum and double integrate it.

"""
# %%
# Load EPR Data
# -------------
# First, load DNPLab,

import dnplab as dnp

# %%
# and then import an EPR spectrum. DNPLab can handle spectra recorded on different spectrometers such as the Bruker ElexSys or EMX system. 
data = dnp.dnpImport.load("../data/bes3t/1D_CW.DTA")

# %%
# Next, create a workspace and add the loaded data
ws = dnp.create_workspace()
ws.add("raw", data)

# %%
# Copy the data to the processing buffer 'proc'
ws.copy("raw", "proc")

# %%
# Raw data are now in both ws["proc"] and ws["raw"]. The following processing steps will modify ws["proc"] but leave ws["raw"] untouched so that you can always return
# to the original data. To see what attributes are loaded with the data use:
print("BEFORE PROCESSING: " + str(ws["raw"].attrs.keys()))

# %%
# Process EPR Data and Plot Results
# ---------------------------------
# In this section, we will demonstrate how to process the EPR data and plot the results. Start with renaming the x-axis to NMR style notation to more seamlessly interact with some functions. "f2" is the default for argument dim.
# ws["proc"].rename("G", "f2")

# %%
# You can display the current axis names with
print("raw axis names: " + str(ws["raw"].dims))
print("proc axis names: " + str(ws["proc"].dims))

# %%
# Now, let's perform a baseline correction using a zeroth order polynomial to remove a DC offset
# print(ws["proc"].dims)

dnp.dnpTools.baseline(ws, dim = "B0", type = "polynomial", order = 0)

# %%
# Keep a copy of the processed spectrum before integrating
ws.copy("proc", "spec")

# %%
# Now let's double integrate the EPR spectrum. 
dnp.dnpTools.integrate(ws, dim = "B0", type = "double")

# %%
# If needed, access your processed spectrum as follows:

# print(ws["proc"].keys())

x_axis = ws["raw"].coords["B0"]                          # Imported field axis
spectrum = ws["spec"].values                            # Spectrum after any processing steps, before integration
# first_integral = ws["proc"].attrs["first_integral"]     # First integral
double_integral = ws["integrals"].values                     # Double integral
print(double_integral)


# %%
# Your originally imported spectrum is still accessible using ws["raw"].values.
#
# To see if any attributes have been added after processing
print("AFTER PROCESSING: " + str(ws["proc"].attrs.keys()))

# %%
# First, let's start with plotting the original EPR spectrum (note the 'raw' keyword) using dnpResults
dnp.dnpResults.plot(ws["raw"])
dnp.dnpResults.plt.title("EPR Spectrum, raw data")
dnp.dnpResults.plt.xlabel("Magnetic Field (mT)")
dnp.dnpResults.plt.grid(True)
dnp.dnpResults.plt.tight_layout()
dnp.dnpResults.show()

# %%
# Note the DC offset of about -0.5.


# %%
# Hint: Use matplotlib to plot arrays such as the attributes "first_integral" or "baseline" over the spectrum, in this case access the spectrum directly as '.values'

# import matplotlib.pyplot as plt


# plt.figure
# plt.plot(ws["raw"].coords["B0"], ws["spec"].values, label = "Corrected EPR Spectrum")
# plt.plot(ws["raw"].coords["B0"], ws["raw"].values, label = "Original EPR Spectrum")
# plt.plot(ws["raw"].coords["B0"], ws["proc"].attrs["baseline"], label = "Baseline")
# plt.plot(ws["raw"].coords["B0"], ws["integrals"].values, label = "DI")
# # plt.plot(
# #     ws["raw"].coords["B0"], ws["proc"].attrs["first_integral"], label = "First Integral"
# # )
# plt.legend()
# plt.title("EPR Spectrum, Baseline Correction, Double Integration")
# plt.grid(True)
# plt.tight_layout()
# plt.show()
