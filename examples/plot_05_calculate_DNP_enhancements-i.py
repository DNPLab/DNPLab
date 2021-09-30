# %%
"""
.. _05-calculate-dnp-enhancements-i:

===================================
05 - Calculate DNP Enhancements (I)
===================================

This example demonstrates how to import DNP-NMR data from an hdf5 file, calculate the dnp enhancements and plot it. If you don't know how to create the workspace yet, please take a look at this tutorial first: :ref:`04-create-a-2d-dnpdata-object-from-individual-spectra`.
"""
# %%
# Load NMR Spectra
# ----------------
# In this example, we will calculate the enhancement for each DNP spectrum, and create a figure of the DNP enhancement versus the microwave power. We will import the 2D dnpdata object created in the previous sample. If you are not yet familiar with how to concatenate individual spectra into the 2D danpdata object, check out this tutorial: :ref:`04-create-a-2d-dnpdata-object-from-individual-spectra`
#
# First, load the 2D dnplab data object and assign it to a workspace (here ws). Once the data is loaded the workspace will have 2 dnpdata objects, the raw data ("raw") and the processed NMR spectra ("proc").

import dnplab as dnp

file_name_path = "../data/h5/PowerBuildUp.h5"
ws = dnp.dnpImport.load(file_name_path)

# %%
# Calculate DNP Enhancement Factors
# ---------------------------------
# The DNP enhancement factors are calculated using the "calculate_enhancement" function. 
dnp.dnpTools.integrate(ws)

# %%% Write note about integrate


dnp.dnpNMR.calculate_enhancement(ws)

# %%
# .. note::
#     The default behavior of the ``calculate_enhancement`` function is to use the first spectrum as the Off signal. If this is the case, the argument ``off_spectrum`` is not necessary unless you want to specify the slice that contains the off spectrum.
#     The ``calculate_enhancement``` function can also calculate the enhancement for specific regions of the spectrum. THis behavior will be discussed in the next example (:ref:`07_align_nmr_spectra`).

# %%
# If needed, access your array of enhancements as:
enhancements = ws["enhancements"].values

# %%
# Plot Enhancement Data
# ---------------------
# Finally, we can plot the enhancement data versus the microwave power.

dnp.dnpResults.figure()
dnp.dnpResults.plot(ws["enhancements"],linestyle = '-', marker = 'o', fillstyle = 'none')
dnp.dnpResults.plt.xlabel("Microwave Power (dBm)")
dnp.dnpResults.plt.ylabel("ODNP Enhancement Factor")
dnp.dnpResults.plt.title("10 mM TEMPO in Toluene")
dnp.dnpResults.plt.grid(True)
dnp.dnpResults.show()
