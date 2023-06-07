# %%
"""
.. _plot_01_align_nmr_spectra:

=================
Align NMR Spectra
=================

In example :ref:`plot_01_load_2D_calculate_DNP_enhancements` we calculated the ODNP enhancement for all protons in toluene, across the entire spectrum. However, the resolution of the spectrum is high enough to resolve the individual NMR peaks of the methyl group and the aromatic protons. This will allow us to calculate the enhancements for the individual peaks. However, as common in low-field NMR spectroscopy, we first have to correct spectra because of the field drift. While magnetic field drift is not an issue in high-field NMR systems using a superconducting magnet, this is not the case in low field systems that either use a electromagnet or a permanent magnet. Peaks will drift and simply averaging over a longer period of time will result in broadened peaks and therefore decreased resolution.

The align routine implemented in DNPLab overlays each individual scan by maximizing the FFT cross-correlation of each scan with a reference scan. This method shifts the spectra in discrete points and calculates the cross-correlation function between the two spectra for each step. The maximum correlation corresponds to the optimum overlap (maximum of the cross-correlation function) between the spectra. This method has successfully been used in `ODNP Spectroscopy <https://linkinghub.elsevier.com/retrieve/pii/S1090780720300379>`_. For more information about alignment methods for NMR spectra check out the article by `Vu et al. <https://doi.org/10.3390/metabo3020259>`_.

"""
# %%
# Load NMR Spectral Data
# ----------------------
# Start with importing the NMR spectra from the previous example (:ref:`plot_04_create_dnpdata_object_from_individual_files`).

import dnplab as dnp

file_name_path = "../../data/h5/PowerBuildUp.h5"
data = dnp.load(file_name_path)
data.attrs["experiment_type"] = "nmr_spectrum"

# %%
# When plotting the data notice that the peaks are not aligned due, in this case, to magnetic field drift.

dnp.fancy_plot(data, xlim=[-10, 20], title="ODNP, 10 mM TEMPO in Toluene, Not Aligned")
dnp.plt.show()

# %%
# Align Spectra
# -------------
# To align the NMR spectra call the align function ``dnp.align()``. This process may take a couple of seconds if the 2D NMR data set is large.

data_aligned = dnp.ndalign(data)

# # %%
# # Next, let's plot the aligned data. Note, that the spectrum is not referenced correctly.

dnp.fancy_plot(
    data_aligned, xlim=[-10, 20], title="ODNP, 10 mM TEMPO in Toluene, Aligned"
)
dnp.plt.show()

# %%
# Save Aligned Spectra
# --------------------
# The last step is to save the aligned spectra for further analysis. Similar to the previous example (:ref:`plot_04_create_dnpdata_object_from_individual_files`) we will save the entire workspace in the h5 file format.

file_name_path = "../../data/h5/PowerBuildUpAligned.h5"
dnp.save(data_aligned, file_name_path, overwrite=True)
