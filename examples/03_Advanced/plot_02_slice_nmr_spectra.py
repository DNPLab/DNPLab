# %%
"""
.. _plot_02_slice_nmr_spectra:

=================
Slice NMR Spectra
=================

This example demonstrates how to slice or select range of NMR spectra. It might be helpful to further data processing.

"""
# %%
# Load NMR Spectral Data
# ----------------------
# Start with importing and plotting the NMR spectra from the previous example (:ref:`plot_04_create_dnpdata_object_from_individual_files`).

import dnplab as dnp

file_name_path = "../../data/h5/PowerBuildUp.h5"
data = dnp.load(file_name_path)
data.attrs["experiment_type"] = "nmr_spectrum"
dnp.fancy_plot(data, xlim=[-10, 20], title="ODNP, 10 mM TEMPO in Toluene")
dnp.plt.show()

# %%
# Slice Spectra Using Integers
# -------------
# In some situations, it is useful to select one or multiple spectra to plot data or to perform further data processing. 
# We will select the a range of data from index 0 to 4. Indexing works similar to other Python objects. However, you must specify the name of the dimension which you want to index.

multiple_spectra_data = data['Power', 0:5] 

# # %%
# # Next, let's plot the sliced spectrum. Note, that the spectrum is not referenced correctly.

dnp.fancy_plot(
    multiple_spectra_data, xlim=[-10, 20], title="ODNP, 10 mM TEMPO in Toluene, Slice Spectrum Using Integers"
)
dnp.plt.show()

# %%
# Slice Spectra Using Floats
# -------------
# In other situations, using the float numbers in index will result in the selection in coordinates.
# In this example, we slice the data range to -8 and 18 ppm.

multiple_spectra_data = data['f2', (-8., 18.)]

# # %%
# # Next, let's plot the sliced data. Note, that the spectrum is not referenced correctly.

dnp.fancy_plot(
    multiple_spectra_data, xlim=[-10, 20], title="ODNP, 10 mM TEMPO in Toluene, Slice Spectra Using Floats"
)
dnp.plt.show()

# %%
# Slice Spectra with Selection of Multiple Dimension
# -------------
# If the data has more than 2 dimension and you want to select dataset, you can specify multiple names of the dimension and the range you want to index. Eg. data['x', 1:10, 'y', :, 'z', (3.5, 7.5)].
# In this example, we select only one spectrum with range between -8. and 18.

multiple_dimension_data = data['f2', (-8., 18.), 'Power', 0]

# # %%
# # Next, let's plot the sliced data. Note, that the spectrum is not referenced correctly.

dnp.fancy_plot(
    multiple_dimension_data, xlim=[-10, 20], title="ODNP, 10 mM TEMPO in Toluene, One Spectra between -8. and 18. left"
)
dnp.plt.show()