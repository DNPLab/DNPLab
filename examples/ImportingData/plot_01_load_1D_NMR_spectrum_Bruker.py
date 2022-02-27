# %%
"""
.. _01_load_1d_spectrum_topspin:

======================================
Load 1D NMR spectrum in TopSpin format
======================================

This example demonstrates how to load and process a single 1D NMR spectrum recorded on a Bruker spectrometer acquired using TopSpin. This examples also demonstrates how to access the various components of the data object and generally how to interact with the functions.

"""
# %%
# First you need to prepare the Python environment by importing DNPLab. In this example we also have to imp, also numpy is imported. However, this is only necessary if you would like to access numpy functions and procedures directly.

import dnplab as dnp
import matplotlib.pylab as plt

# %%
# Load and Process the NMR Spectrum
# ---------------------------------
# In this example we use a 1D NMR spectrum aqcuired using TopSpin. Example data is located in the data folder. All data enters into the same object structure so this example applies to any NMR format. DNPLab will try to identify the format of the NMR spectrum (e.g. TopSpin, VnmrJ, Kea, etc.). This will work in most cases. If the autodetect fails, the format can be explicitely given using the data_type attribute of the load function (e.g. data_type = "topspin").

data = dnp.load("../../data/topspin/1")

# %%
# The spectral data is now stored in the data object. DNPLab uses object oriented programming, therefore, data is not just a variable. It is a dnpData object. Next, we will perform some basic processing, methods that are pretty standard in NMR spectroscopy. For this just pass the dnpData object from one function to the next.

# data = dnp.remove_background(data)
data = dnp.apodize(data, lw = 100)
data = dnp.fourier_transform(data)

# %%
# This procedure removes the DC offset, performs an apodization with a window function (100 Hz linewidth), followed by a Fourier transforms. By default a zero filling of factor two is performed prioer to Fourier transforming the spectrum.
#
# NMR spectra are typically phase sensitive (complex data) and the final processing step is to phase the spectrum. In DNPLab you can either manually phase correct the spectrum using the 'autophase' function using the 'manual' keyword and by giving a phase explicitely:

data = dnp.autophase(data, method = "manual", phase = -1*dnp.constants.pi / 2)

# %%
# Or by using the autophase function to automatically search for and apply the best phase angle

data = dnp.autophase(data, force_positive = True)

# %%
# Plot Processed Data
# -------------------

# %%
# DNPLab automatically takes care of the name of the coordinates. For NMR spectra, the dimensions are named "t2" and "t1" for the direct and indirect dimension of a FID, respectively. After the Fourier Transformation these dimensions will be renamed to "f2" and "f1", respectively. Dimensions can be renamed by the user. However, several features in DNPLab take advantage of the dimension name and processing or plotting functions will change their behavior according to the type of the spectrum.

# %%
# To plot the NMR spectrum simply use the following command. By default, only the real part of the NMR spectrum is shown.

dnp.dnplabplot(data)

# %%
# Note, that the DNPLab plotting function will automatically reverse the direction of the x axis, as it is custom for NRM spectra and axis labels are automatically generated

# %%
# The integrated plotting function allows for additional input arguments. For example to only show a portion of the spectrum use the xlim argument. A title can be passed to the function using the title argument. This can also be a variable. For example:

sampleTag = 'ODNP Experiment of 10 mM TEMPO in Water'

dnp.dnplabplot(data, xlim = [-100, 100], title = sampleTag)

# %%
# Advanced Usage
# --------------
# The dnpData objects stores many more attributes for the data set, which are loaded automatically when loading the experimental data. For example, the NMR frequency can be easily accessed as follows:

nmr_frequency = data.attrs["nmr_frequency"]
print(nmr_frequency)

# %%
# All variables in DNPLab are stored in SI units. Therefore, the NMR frequency is given in (Hz).
# If needed, access your processed spectrum as follows:

x_axis = data.coords["f2"]        # ppm axis
spectrum = data.values            # spectrum