# %%
"""
.. _07_align_nmr_spectra:

======================
07 - Align NMR Spectra
======================

This example demonstrates how to align NMR spectra. While magnetic field drift is not an issue in high-field NMR systems using a superconducting magnet, this is not the case in low field systems that either use a electromagnet or a permanent magnet.

The align routine implemented in DNPLab overlays each individual scan by maximizing the FFT cross-correlation of each scan with a reference scan. This method shifts the spectra in discrete points and calculates the cross-correlation function between the two spectra for each step. The maximum correlation corresponds to the optimum overlap (maximum of the cross-correlation function) between the spectra. This method has successfully been used in `ODNP Spectroscopy <https://linkinghub.elsevier.com/retrieve/pii/S1090780720300379>`_. For more information about alignment methods for NMR spectra check out the article by `Vu et al. <https://doi.org/10.3390/metabo3020259>`_.

"""
# %%
# Load NMR Spectral Data
# ----------------------
# Start with importing data and create the DNPLab workspace
import dnplab as dnp
file_name_path = "../data/topspin/304"
data = dnp.dnpImport.load(file_name_path)

ws = dnp.create_workspace()
ws.add("raw", data)
ws.copy("raw", "proc")

# %%
# Align Spectra
# -------------
# The spectra shown here are the individual NMR spectra for different recovery times. Note, how this spectra are not perfectly aligned. This is can changed using the alining function. A more detailed description of the alining routine of DNPLab is given in Tutorial XXXX.

