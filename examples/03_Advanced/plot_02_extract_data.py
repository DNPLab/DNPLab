# %%
"""
.. _plot_02_extract_data:

==========================
Extract Individual Spectra
==========================

This example demonstrates how to extract data from a multidimensional dnpdata object. For example, extracting a single spectrum (or a set of spectra) for a 2D set of spectra.

"""
# %%
# Load 2D NMR Data
# ----------------
# For this example we will use a two-dimensional ODNP-enhanced NMR data set. Start by importing and plotting the 2D NMR data set created in the previous example (:ref:`plot_01_align_nmr_spectra`).

import dnplab as dnp

file_name_path = "../../data/h5/PowerBuildUpAligned.h5"
data = dnp.load(file_name_path)
data.attrs["experiment_type"] = "nmr_spectrum"
dnp.fancy_plot(data, xlim=[-10, 20], title="ODNP, 10 mM TEMPO in Toluene")
dnp.plt.show()

# %%
# The example data set consists of ODND-enhanced NMR spectra of a sample of 10 mM TEMPO in toluene. In the second dimension the microwave power is varied.

# %%
# Indexing Using Integers
# -----------------------
# Let's start with the simplest way to extract (index) a single spectrum from a 2D NMR data set. The data set imported here are ODNP-enhanced 1H NMR spectra of toluene. One dimension contains the NMR spectrum (``dim = 'f2'``) and the microwave power is increased throughout the second dimension (``dim = 'power'``). To extract data from the dnpdata object, the user has to specify the name of the dimension and the index of the spectrum. This index has to be an integer and follows the Python convention for indexing (0 is the first index, -1 is the last index). For example, the 4th NMR spectrum at a given microwave power can be extracted from the 2D data set by:

single_spectrum = data['Power', 3]

# %%
# Next, we can plot the spectrum (Note, that the spectrum is not referenced correctly).

dnp.fancy_plot(single_spectrum, xlim=[-10, 20])
dnp.plt.show()

# %%
# The above command creates a new dnpdata object named ``single_spectrum``. If you just want to plot a a specific trace, this can be done by indexing the dnpdata object directly:

dnp.fancy_plot(data['Power', -1], xlim=[-10, 20])
dnp.plt.show()

# %%
# Here, the last spectrum (python index -1) of the dnpdata object ``data`` is plotted.

# %%
# In some situations it is useful to select a subset of spectra from the data set e.g. for further data processing. For example, the first 5 spectra of the dnpdata object ``data`` can be selected by:

sub_data = data['Power', 0:5] 

# %%
# Indexing Using Floats
# ---------------------
# In the previous example, a single spectrum or a sub-set of spectra was selected using an integer number to index the data set. This is convenient if this dimension has only a small number of indexes. However, when this dimension is large it is less convenient to use an integer. For example, it is not very convenient to index a 1D NMR spectrum with 8,192 points by an integer indicating the position in the vector/array. Here, it is much more convenient to give the axis value to extract the data that falls within these boundaries.

# %% For example, the NMR spectrum for a chemical shift range from -8 to 18 ppm can be extracted by executing the following command:

sub_data = data['f2', (-8., 18.)]

# %%
# The spectra are plotted below. Note the range of the data set is -8 to 18 ppm. The new dnpdata object ``sub_data`` only contains the selected data region, the remaining data points are discarded.

dnp.fancy_plot(sub_data, xlim=[-15, 25], title="ODNP, 10 mM TEMPO in Toluene, Slice Spectra Using Floats")
dnp.plt.show()

# %%
# Indexing Multiple Dimensions
# ----------------------------
# For multi-dimensional data sets, the user can specify multiple dimensions at once. For example ``data['x', 1:10, 'y', :, 'z', (3.5, 7.5)]``. In the following example we will select only one single NMR spectrum (e.g. ``'Power' , 0``) within a spectral range of -8 to 18 ppm.

sub_data = data['f2', (-8., 18.), 'Power', 0]

dnp.fancy_plot(sub_data)
dnp.plt.show()