# %%
"""
04 - Calculate DNP Enhancements
===============================

This example demonstrates how to import DNP-NMR data and calculate the dnp enhancement.
"""
# %%
# Load NMR Spectra
# ----------------
# For this example a set of 1D NMR spectra is imported, which was recorded using different microwave powers. Since this series was not saved as a 2D spectrum/aray, we first have to create a 2D data object, which will then be processed. DNPLab this task can be performed very efficiently in a single line of code.

import dnplab as dnp
import numpy as np

powers = np.array([1,2,3,4])

filenames = [
    "../data/prospa/toluene_10mM_Tempone/1/data.1d",
    "../data/prospa/toluene_10mM_Tempone/29/data.1d",
    "../data/prospa/toluene_10mM_Tempone/30/data.1d",
    "../data/prospa/toluene_10mM_Tempone/31/data.1d",
]

data = dnp.dnpImport.load(filenames, dim = "power", coord = powers)

# %%
# In the code section shown above, first an array of power values is created. The paths to the individual spectra are stored in the Python list called "filenames". To load the data, the dnpImport.load function used. A new dimension is automatically created called "power" and the data is concatenated into a single 2D dnpData object.
#
# Finally, we can create the workspace, add the data to the "raw" object, and copy the "raw" data to the processing buffer.

ws = dnp.create_workspace()
ws.add("raw", data)
ws.copy("raw", "proc")

# %%
# Porcess the NMR Spectra and Calculate Enhancement Factors
# ---------------------------------------------------------
# Once the 2D data set is created, NMR processing is straightforward.
dnp.dnpNMR.remove_offset(ws)
dnp.dnpNMR.window(ws, linewidth = 10)
dnp.dnpNMR.fourier_transform(ws, zero_fill_factor = 2)

# The DNP enhancement factors can be calculated. 
dnp.dnpNMR.calculate_enhancement(ws, off_spectrum = 0)

# # %% [markdown]
# # If needed, access your array of enhancements as:
# enhancements = ws["enhancement"].values
# # %%

# %%
# Note: The default behavior of the calculate_enhancement function is to use the first spectrum as the Off singal and the argument off_spectrum is not necessary unless the off spectrum is not the first slice. 






# %%
# Fit Enhancement Factors
# -----------------------
# 
# 

# %%
# Plot Enhancement Data and Fit
# -----------------------------
# 
# 

dnp.dnpResults.figure()
dnp.dnpResults.plot(ws["enhancement"])
# dnp.dnpResults.legend()
dnp.dnpResults.plt.grid(True)
dnp.dnpResults.show()

