# %%
"""
.. _07_align_nmr_spectra:

======================
07 - Align NMR Spectra
======================

In example :ref:`05-calculate-dnp-enhancements-i` we calculated the enhancement for the toluene protons across the entire spectrum. However, in this particular case, the resolution is high enough to resolve the NMR peaks of the methyl group and the aromatic protons. To calculate the enhancements for the individual peaks first the spectra need to be aligned.

While magnetic field drift is not an issue in high-field NMR systems using a superconducting magnet, this is not the case in low field systems that either use a electromagnet or a permanent magnet. Peaks will drift and simply averaging over a longer period of time will result in broadened peaks and therefore decreased resolution.

The align routine implemented in DNPLab overlays each individual scan by maximizing the FFT cross-correlation of each scan with a reference scan. This method shifts the spectra in discrete points and calculates the cross-correlation function between the two spectra for each step. The maximum correlation corresponds to the optimum overlap (maximum of the cross-correlation function) between the spectra. This method has successfully been used in `ODNP Spectroscopy <https://linkinghub.elsevier.com/retrieve/pii/S1090780720300379>`_. For more information about alignment methods for NMR spectra check out the article by `Vu et al. <https://doi.org/10.3390/metabo3020259>`_.

"""
# %%
# Load NMR Spectral Data
# ----------------------
# Start with importing the NMR spectra from the prevous example (:ref:`04-create-a-2d-dnpdata-object-from-individual-spectra`). 

import dnplab as dnp

file_name_path = "../data/h5/PowerBuildUp.h5"
ws = dnp.dnpImport.load(file_name_path)

# %%
# When plotting the data you notice that the peaks are not on top of each other due to the drift of the magnetic field.

dnp.dnpResults.figure()
dnp.dnpResults.plot(ws["proc"].real)
dnp.dnpResults.xlim([30, -30])
dnp.dnpResults.plt.xlabel("Chemical Shift (ppm)")
dnp.dnpResults.plt.ylabel("Signal Amplitude (a.u.)")
dnp.dnpResults.plt.title("DNP-NMR Spectra of 10 mM TEMPO in Toluene, Not Aligned")
dnp.dnpResults.plt.grid(True)
dnp.dnpResults.show()

# %%
# Align Spectra
# -------------
# The spectra shown in the figure above are not aligned. Next, we will use the align function to correct for the magnetic field drift. First we start with making a copy of the processing buffer ("proc") and saving it under the name "NotAligned". In the next step the data is aligned by calling the function ``dnp.dnpNMR.align()``. This process will take a couple of seconds, especially if the 2D NMR data set is large. Once this operation is completed, we copy the data set from the processing buffer to a new dnpdata object simply called "Aligned".

ws.copy("proc", "NotAligned")                       # Copy processing buffer
dnp.dnpNMR.align(ws)                                # Align NMR spectra
ws.copy("proc", "Aligned")                          # Copy aligned spectra to new object

# %%
# Next, let's plot the aligned data.

ws["Aligned"].coords["f2"] += 4.4

dnp.dnpResults.figure()
dnp.dnpResults.plot(ws["Aligned"].real)
dnp.dnpResults.xlim([20, 0])
dnp.dnpResults.plt.xlabel("Chemical Shift (ppm)")
dnp.dnpResults.plt.ylabel("Signal Amplitude (a.u.)")
dnp.dnpResults.plt.title("DNP-NMR Spectra of 10 mM TEMPO in Toluene, Aligned")
dnp.dnpResults.plt.grid(True)
dnp.dnpResults.show()

# %%
# Save Aligned Spectra
# --------------------
# The last step is to save the aligned spectra. Similar to the previous example (:ref:`04-create-a-2d-dnpdata-object-from-individual-spectra`) we will save the entire workspace in the hdf5 file format.

file_name_path = "../data/h5/PowerBuildUpAligned.h5"
dnp.dnpSave.save(ws,file_name_path, overwrite = True)
