# %%
"""
04 - Calculate DNP Enhancements
===============================

This example demonstrates how to import DNP-NMR data and calculate the dnp enhancement.
"""
# %%
# Load NMR Spectra
# ----------------
# For this example a set of 1D NMR spectra is imported, which was recorded using different microwave powers. Since this series was not saved as a 2D spectrum, we first have to create a 2D data object, which will then be processed.

import dnplab as dnp
import numpy as np

temp_data_list = []
powers = np.array([1,2,3,4])

filenames = [
    "../data/prospa/toluene_10mM_Tempone/1/data.1d",
    "../data/prospa/toluene_10mM_Tempone/29/data.1d",
    "../data/prospa/toluene_10mM_Tempone/30/data.1d",
    "../data/prospa/toluene_10mM_Tempone/31/data.1d",
]

for filename in filenames:
    tmp = dnp.dnpImport.load(filename, data_type = "prospa")
    temp_data_list.append(tmp)

data = dnp.concat(temp_data_list, "power", powers)

# %%
# In the code section shown above, first an empty list (temp_data_list) is created and an array of power values. The paths to the individual spectra is stored in the Python list called filenames. To load the data, we loop over each entry of the Python list and add the loaded data to a python list. To create the DNPLab data object, the list is then concatenated.
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


dnp.dnpNMR.calculate_enhancement(ws)




print(ws["enhancement"])


dnp.dnpResults.figure()
# dnp.dnpResults.plot(ws["enhancement"])
dnp.dnpResults.plot(ws["proc"])

# dnp.dnpResults.legend()
# dnp.dnpResults.plt.grid(True)
dnp.dnpResults.show()



    # enhancements.append(ws_on["enhancement"].values)

# enhancements = []
    # ws.add("proc", data)  # add directly to the processing buffer "proc"



