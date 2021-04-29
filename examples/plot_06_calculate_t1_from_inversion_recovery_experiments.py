# %%
"""
.. _06-calculate_t1_from_inversion_recovery:

=====================================================
06 - Calculate T1 from Inversion-Recovery Experiments
=====================================================

This example demonstrates how to import data from an inversion recovery NMR experiment and determine the T1 relaxation rate through a fit.
"""
# %%
# Load Inversion Recovery Spectra
# -------------------------------
# Start with importing data and create the DNPLab workspace
import dnplab as dnp
file_name_path = "../data/topspin/304"
data = dnp.dnpImport.load(file_name_path)

ws = dnp.create_workspace()
ws.add("raw", data)
ws.copy("raw", "proc")

# %%
# Next, we process the FIDs to obtain the NMR spectrum
dnp.dnpNMR.remove_offset(ws)
dnp.dnpNMR.window(ws, linewidth = 10)
dnp.dnpNMR.fourier_transform(ws, zero_fill_factor = 2)

dnp.dnpNMR.autophase(ws)

dnp.dnpResults.plot(ws["proc"].real)
dnp.dnpResults.xlim([-50, 30])
dnp.dnpResults.plt.xlabel("Chemical Shift [ppm]")
dnp.dnpResults.plt.ylabel("Signal Amplitude [a.u.]")
dnp.dnpResults.plt.title("Inversion Recovery")
dnp.dnpResults.plt.grid(True)
dnp.dnpResults.show()

# %%
# Align Inversion Recovery Spectra
# --------------------------------
# The spectra shown here are the individual NMR spectra for different recovery times. Note, how this spectra are not perfectly aligned. This is can changed using the alining function. A more detailed description of the alining routine of DNPLab is given in tutorial :ref:`07_align_nmr_spectra`.

dnp.dnpNMR.align(ws)

dnp.dnpResults.plot(ws["proc"].real)
dnp.dnpResults.xlim([-50, 30])
dnp.dnpResults.plt.xlabel("Chemical Shift [ppm]")
dnp.dnpResults.plt.ylabel("Signal Amplitude [a.u.]")
dnp.dnpResults.plt.title("Inversion Recovery, Spectra Aligned")
dnp.dnpResults.plt.grid(True)
dnp.dnpResults.show()

# %%
# To determine the T1 relaxation time, we first integrate the peak intensity across the entire spectrum. To save the spectral data, we move the processed NMR spectrum within the workspace out of the processing buffer into a new buffer called "ft". If we don't do that, the processing buffer (proc) will be overwritten by the integrate function.

ws.copy("proc", "ft")
dnp.dnpTools.integrate(ws)

# %%
# Display Results
# ---------------

# # Fit inversion-recovery build-up, display T1 and plot results

dnp.dnpFit.exponential_fit(ws, type="T1")
print("T1 value (sec) = " + str(ws["fit"].attrs["T1"]))

dnp.dnpResults.figure()
dnp.dnpResults.plot(ws["integrals"], "o", fillstyle="none")
dnp.dnpResults.plot(ws["fit"])
dnp.dnpResults.plt.xlabel("Time t1 [s]")
dnp.dnpResults.plt.ylabel("Signal Amplitude [a.u.]")
dnp.dnpResults.plt.title("Inversion Recovery")
dnp.dnpResults.plt.grid(True)
dnp.dnpResults.show()
