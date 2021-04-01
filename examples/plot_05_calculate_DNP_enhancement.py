# %%
"""
.. _05-calculate-dnp-enhancements:

===============================
05 - Calculate DNP Enhancements
===============================

This example demonstrates how to import DNP-NMR data and calculate the dnp enhancement.
"""
# %%
# Load NMR Spectra
# ----------------
# In this example, we will calculate the enhancement for each DNP spectrum, and create a figure of the DNP enhancement versus the microwave power. We will import the 2D dnpdata object created in the previous sample. If you are not yet familiar with how to concatenate individual spectra into the 2D danpdata object, check out this tutorial: :ref:`04-create-a-2d-dnpdata-object-from-individual-spectra`
#
# First, load the 2D data object. This will load the workspace (ws) which will have the raw data ("raw") and the processed NMR spectra ("proc").

import dnplab as dnp
import numpy as np

# file_name_path = "../data/prospa/toluene_10mM_Tempone/PowerBuildUp.h5"
# dnp.dnpImport.load(file_name_path)

powers = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40])

filenames = [
    "../data/prospa/toluene_10mM_Tempone/1/data.1d",
    "../data/prospa/toluene_10mM_Tempone/2/data.1d",
    "../data/prospa/toluene_10mM_Tempone/3/data.1d",
    "../data/prospa/toluene_10mM_Tempone/4/data.1d",
    "../data/prospa/toluene_10mM_Tempone/5/data.1d",
    "../data/prospa/toluene_10mM_Tempone/6/data.1d",
    "../data/prospa/toluene_10mM_Tempone/7/data.1d",
    "../data/prospa/toluene_10mM_Tempone/8/data.1d",
    "../data/prospa/toluene_10mM_Tempone/9/data.1d",
    "../data/prospa/toluene_10mM_Tempone/10/data.1d",
    "../data/prospa/toluene_10mM_Tempone/11/data.1d",
    "../data/prospa/toluene_10mM_Tempone/12/data.1d",
    "../data/prospa/toluene_10mM_Tempone/13/data.1d",
    "../data/prospa/toluene_10mM_Tempone/14/data.1d",
    "../data/prospa/toluene_10mM_Tempone/15/data.1d",
    "../data/prospa/toluene_10mM_Tempone/16/data.1d",
    "../data/prospa/toluene_10mM_Tempone/17/data.1d",
    "../data/prospa/toluene_10mM_Tempone/18/data.1d",
    "../data/prospa/toluene_10mM_Tempone/19/data.1d",
    "../data/prospa/toluene_10mM_Tempone/20/data.1d",
    "../data/prospa/toluene_10mM_Tempone/21/data.1d",
    "../data/prospa/toluene_10mM_Tempone/22/data.1d",
    "../data/prospa/toluene_10mM_Tempone/23/data.1d",
    "../data/prospa/toluene_10mM_Tempone/24/data.1d",
    "../data/prospa/toluene_10mM_Tempone/25/data.1d",
    "../data/prospa/toluene_10mM_Tempone/26/data.1d",
    "../data/prospa/toluene_10mM_Tempone/27/data.1d",
    "../data/prospa/toluene_10mM_Tempone/28/data.1d",
    "../data/prospa/toluene_10mM_Tempone/29/data.1d",
    "../data/prospa/toluene_10mM_Tempone/30/data.1d",
    "../data/prospa/toluene_10mM_Tempone/31/data.1d",
    "../data/prospa/toluene_10mM_Tempone/32/data.1d",
    "../data/prospa/toluene_10mM_Tempone/33/data.1d",
    "../data/prospa/toluene_10mM_Tempone/34/data.1d",
    "../data/prospa/toluene_10mM_Tempone/35/data.1d",
    "../data/prospa/toluene_10mM_Tempone/36/data.1d",
    "../data/prospa/toluene_10mM_Tempone/37/data.1d",
    "../data/prospa/toluene_10mM_Tempone/38/data.1d",
    "../data/prospa/toluene_10mM_Tempone/39/data.1d",
    "../data/prospa/toluene_10mM_Tempone/40/data.1d",
    "../data/prospa/toluene_10mM_Tempone/41/data.1d",
]

data = dnp.dnpImport.load(filenames, dim = "power", coord = powers)
ws = dnp.create_workspace()
ws.add("raw", data)
ws.copy("raw", "proc")

dnp.dnpNMR.remove_offset(ws)
dnp.dnpNMR.window(ws, linewidth = 10)
dnp.dnpNMR.fourier_transform(ws, zero_fill_factor = 2)

# %%
# Calculate DNP Enhancement Factors
# ---------------------------------
# The DNP enhancement factors are calculated using the "calculate_enhancement" function. 
dnp.dnpNMR.calculate_enhancement(ws, off_spectrum = 0)

# %% [markdown]
# If needed, access your array of enhancements as:
enhancements = ws["enhancement"].values
# %%
# Note: The default behavior of the calculate_enhancement function is to use the first spectrum as the Off singal and the argument "off_spectrum" is not necessary unless the off spectrum is not the first slice.

# %%
# Plot Enhancement Data and Fit
# -----------------------------
# Finally, we can plot the enhancement data and the corresponding fit. 

dnp.dnpResults.figure()
dnp.dnpResults.plot(ws["enhancement"],linestyle = '-', marker = 'o', fillstyle = 'none')
dnp.dnpResults.plt.xlabel("Microwave Power (dBm)")
dnp.dnpResults.plt.ylabel("ODNP Enhancement Factor")
dnp.dnpResults.plt.title("10 mM TEMPO in Toluene")
dnp.dnpResults.plt.grid(True)
dnp.dnpResults.show()
