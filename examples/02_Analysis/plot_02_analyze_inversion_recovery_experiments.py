# %%
"""
.. _plot_02_analyse_inversion_recovery_experiments:

=========================================
Analyze T1 Inversion-Recovery Experiments
=========================================

This example demonstrates how to import TopSpin data from an inversion recovery NMR experiment and determine the T1 relaxation rate through a fit.
"""
# %%
# Load Inversion Recovery Spectra
# -------------------------------
# Start with importing data and creating the DNPLab workspace. Note that the data are added to the workspace as 'raw' then copied to the processing buffer 'proc' to preserve the raw data during processing.
import dnplab as dnp

file_name_path = "../../data/topspin/304"
data = dnp.load(file_name_path)
data.attrs["experiment_type"] = "nmr_spectrum"

# %%
# Next, we process the FIDs to obtain the frequency domain NMR spectrum

# data = dnp.remove_background(data)

# dnp.remove_background(data)
data = dnp.apodize(data, lw=100)
data = dnp.fourier_transform(data)

data = dnp.autophase(data)

# dnp.dnplabplot(data, xlim=[-50, 30], title="Inversion Recorvery")
# dnp.plt.show()

# %%
# Align Inversion Recovery Spectra
# --------------------------------
# The spectra shown here are the individual NMR spectra for different recovery times. Note how the spectra are not perfectly aligned. This is can be fixed using the aligning function. A more detailed description of the aligning routine of DNPLab is given in the tutorial :ref:`07_align_nmr_spectra`.

data = dnp.align(data)

# dnp.dnplabplot(data, xlim=[-50, 30], title="Inversion Recovery, aligned")
# dnp.plt.show()


# %%
# Integrate Spectra
# -----------------
# To determine the T1 relaxation time, we first integrate the peak intensity across the entire spectrum. After integration the workspace will have a new dnpdata object called "integrals" where the integral values and indirect axis are kept.

integrals = dnp.integrate(data)
integrals.attrs["experiment_type"] = "inversion_recovery"

dnp.fancy_plot(integrals)
dnp.plt.show()

# %%
# Fit Inversion Recovery Buildup
# ------------------------------
# Fit the inversion recovery build-up and display the T1 value. After the fit the workspace will have a new dnpdata object called "fit" where the fit array and fit axis are kept.

# dnp.dnpFit.exponential_fit(ws, type="T1")
# print("T1 value (sec) = " + str(ws["fit"].attrs["T1"]))

# dnp.dnpResults.figure()
# dnp.dnpResults.plot(ws["integrals"], "o", fillstyle="none")
# dnp.dnpResults.plot(ws["fit"])
# dnp.dnpResults.plt.xlabel("Time t1 [s]")
# dnp.dnpResults.plt.ylabel("Signal Amplitude [a.u.]")
# dnp.dnpResults.plt.title("Inversion Recovery")
# dnp.dnpResults.plt.grid(True)
# dnp.dnpResults.show()
