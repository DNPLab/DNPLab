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

dnp.fancy_plot(data, xlim=[-50, 30], title="Inversion Recorvery")
dnp.plt.show()

# %%
# Align Inversion Recovery Spectra
# --------------------------------
# The spectra shown here are the individual NMR spectra for different recovery times. Note how the spectra are not perfectly aligned. This is can be fixed using the aligning function. A more detailed description of the aligning routine of DNPLab is given in the tutorial :ref:`plot_01_align_nmr_spectra`.

data = dnp.align(data)

# dnp.dnplabplot(data, xlim=[-50, 30], title="Inversion Recovery, aligned")
# dnp.plt.show()


# %%
# Integrate Spectra
# -----------------
# To determine the T1 relaxation time, we first integrate the peak intensity across the entire spectrum. After integration the workspace will have a new dnpdata object called "integrals" where the integral values and indirect axis are kept. To use ``fancy_plot`` the attribute experiment_type has to be changed to "inversion_recovery".

integrals = dnp.integrate(data)
integrals.attrs["experiment_type"] = "inversion_recovery"

dnp.fancy_plot(integrals)
dnp.plt.show()


# %%
# Fit Inversion-Recovery Experiment
# ---------------------------------
# Once the peaks are integrated the T1 value can be determined by fitting the signal intensities to a function that describes the inversion-recovery experiment.


out = dnp.inversion_recovery_fit(integrals)

# fit = out['fit']
# popt = out['popt']
# err = out['err']

# dnp.fancy_plot(integrals)
# dnp.plot(fit, '-')

# # print(popt)


# print(T1.values)

# dnp.plt.show()
