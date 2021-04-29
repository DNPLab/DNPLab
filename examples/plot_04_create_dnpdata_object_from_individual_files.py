# %%
"""
.. _04-create-a-2d-dnpdata-object-from-individual-spectra:

=======================================================
04 - Create a 2D dnpdata Object from Individual Spectra
=======================================================

This example demonstrates how to import a list of DNP-NMR spectra and create a 2D dnpdata object.

Depending on how you record a set of DNP-NMR experiments, you will either end up with a single file corresponding to a 2D array of spectra (in which case you can skip to the next example ...) or with a set of individual files. A common example is recording the DNP enhancement as a function of the microwave power. For easy data handling, these individual spectra can be concatanated in a single dnpdata object for easy processing and analyzing of the data.
"""
# %%
# Load NMR Spectra
# ----------------
# For this example a set of 1D NMR spectra is imported. Each spectrum is recorded using a different microwave power. The import function of DNPLab can handle a list of spectra and will automatically create the dnpdata object. To load multiple spectra first create a list of paths to the individual spectra (alternatively, you can loop over the folder index, however, for educational purposes we keep this simple for now).

import dnplab as dnp
import numpy as np

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

powers = np.array(
    [
        0,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        28,
        29,
        30,
        31,
        32,
        33,
        34,
        35,
        36,
        37,
        38,
        39,
        40,
    ]
)  # Power in dBm

data = dnp.dnpImport.load(filenames, dim="Power", coord=powers)

# %%
# In the next step create an array with the power levels. The length of this array should match the number of spectra. The Python list "filenames" and the array of power levels will become input arguments to the load function. Here, the dimension is called "Power" and the values stored in "powers" serves as the "coord" input argument. When importing the spectra DNPLab will automatically create a 2D object with a new dimension namend "Power" and the data is concatenated into a single 2D dnpdata object.
#
# Finally, we can create the workspace, add the data to the "raw" object, and copy the "raw" data to the processing buffer.

ws = dnp.create_workspace()
ws.add("raw", data)
ws.copy("raw", "proc")

# %%
# Process and Save the NMR Spectra
# --------------------------------
# Once the 2D data set is created, NMR processing is straightforward. Here, we apply a line-broadening of 10 Hz, perform a Fourier Transformation, and zero-filling of the data set to twice the number of points.
dnp.dnpNMR.remove_offset(ws)
dnp.dnpNMR.window(ws, linewidth=10)
dnp.dnpNMR.fourier_transform(ws, zero_fill_factor=2)

# %%
# Finally, the 1D spectra are plotted.
dnp.dnpResults.figure()
dnp.dnpResults.plot(ws["proc"].real)
dnp.dnpResults.xlim([30, -30])
dnp.dnpResults.plt.xlabel("Chemical Shift [ppm]")
dnp.dnpResults.plt.ylabel("Signal Amplitude [a.u.]")
dnp.dnpResults.plt.title("DNP Enhancement Power Build-Up, 10 mM TEMPO in Toluene")
dnp.dnpResults.plt.grid(True)
dnp.dnpResults.show()

# %%
# Saving the Processed Data
# -------------------------
# DNPLab has built-in capabilities to save large data sets, so we can save the already concatenated and processed NMR data in a single file and load just this file for further processing.

file_name_path = "../data/h5/PowerBuildUp.h5"
dnp.dnpSave.save(ws, file_name_path, overwrite=True)

# %%
# DNPLab saves the 2D dnpdata object in the hdf5 file format. We will use this data in the next example (:ref:`05-calculate-dnp-enhancements-i`) for further processing.
