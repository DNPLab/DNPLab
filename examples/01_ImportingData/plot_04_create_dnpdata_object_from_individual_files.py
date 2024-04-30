# %%
"""
.. _plot_04_create_dnpdata_object_from_individual_files:

=========================================================
Create a 2D dnpdata object from set of individual spectra
=========================================================

This example demonstrates how to import a list of DNP-NMR spectra and create a 2D dnpdata object.

Depending on how you record a set of DNP-NMR experiments, you will either end up with a single file corresponding to a 2D array of spectra (in which case you can skip to the next example ...) or with a set of individual files recorded for example at different microwave power levels. A common example is recording the NMR signal at different levels of microwave power to determine the enhancement at maximum power. Processing each spectrum individually is tedious and time-consuming. To make processing more convenient, the individual NMR spectra can be concatenated in a single dnpdata object for easy processing and analyzing of the data.
"""
# %%
# Load NMR Spectra
# ----------------
# For this example a set of 41 individual 1D NMR spectra are imported. Each spectrum is recorded using a different microwave power. The import function of DNPLab can handle a list of spectra and will automatically create the dnpdata object. To load multiple spectra first create a list of paths to the individual spectra (alternatively, you can loop over the folder index, however, for educational purposes we keep this example simple for now).

import dnplab as dnp
import numpy as np

filepath_prefix = "../../data/prospa/toluene_10mM_Tempone/"

filenames = [filepath_prefix + "%d/data.1d" % i for i in range(1, 42)]

# %%
# In addition, create an array with the power levels. In this example we use numpy to create the array. The length of this array should match the number of spectra. The Python list "filenames" and the array of power levels will become input arguments to the load function. Here, the dimension is called "Power" and the values stored in "powers" serves as the "coord" input argument. When importing the spectra DNPLab will automatically create a 2D object with a new dimension named "Power" and the data is concatenated into a single 2D dnpdata object. In this example the power is given in units of dBm.

powers = np.linspace(0, 40, 41)

# %%
# Now load the data and assign the power array to coord,

data = dnp.load(filenames, dim="Power", coord=powers)

# %%
# Process and Save the NMR Spectra
# --------------------------------
# Once the 2D data set is created, NMR processing is straightforward. Here, we apply a line-broadening of 10 Hz, perform a Fourier Transformation, and zero-filling of the data set to twice the number of points (default of the Fourier transform function).

# dnp.dnpNMR.remove_offset(ws)

data = dnp.apodize(data, lw=10)
data = dnp.fourier_transform(data)

# %%
# Once the raw data are processed it is time to plot the 1D spectra.

sampleTag = "10 mM TEMPO in Toluene"

dnp.plt.figure()
dnp.fancy_plot(data, xlim=[-10, 20])
dnp.plt.title(sampleTag)
dnp.plt.show()
dnp.plt.tight_layout()

# %%
# Saving the Processed Data
# -------------------------
# DNPLab can save large data sets in a single file, so the processed data can be used at a later stage for further processing or analysis.

file_name_path = "../../data/h5/PowerBuildUp.h5"
dnp.save(data, file_name_path, overwrite=True)

# %%
# DNPLab saves the 2D dnpdata object in the hdf5 file format. We will use this data in the next example (:ref:`plot_01_load_2D_calculate_DNP_enhancements`) for further processing.
