# %%
"""
.. _plot_02_load_1D_NMR_spectrum_Kea:

=====================================
Load two 1D NMR spectra in Kea format
=====================================

In this example we demonstrate how to import two ODNP-enhanced NMR spectra, one recorded with a microwave power of 0 W (off-signal) and one with a microwave power of 2 W (on-signal). The spectra are recorded using a Magritek Kea system.

The example script has three different sections:

#. Load and Process Off-Signal
#. Load and Process On-Signal
#. Create a Figure and Plot On/Off Spectra

"""
# %%
# Make sure to start with importing DNPLab
import dnplab as dnp

# %%
# Load and Process Off-Signal
# -----------------------------
# The next section demonstrates how the FID is imported into DNPLab and processed. Processing involves removing any DC offset, followed by a 15 Hz linewidth appodization, prior to performing the Fourier transformation.

########## OFF Signal (P = 0 W) ##########
data_off = dnp.load("../../data/prospa/10mM_TEMPO_Water/1Pulse_20200929/35/data.1d")
data_off.attrs["experiment_type"] = "nmr_spectrum"

data_off = dnp.remove_background(data_off)
data_off = dnp.apodize(data_off, lw=15)
data_off = dnp.fourier_transform(data_off)

# %%
# Load and Process ON-Signal
# ----------------------------
# Importing the on-signal involves the same steps as importing the off-signal. Once processed the data is copied to the results buffer 'onSignal'.

########## ON Signal (P = 2 W) ##########
data_on = dnp.load("../../data/prospa/10mM_TEMPO_Water/1Pulse_20200929/51/data.1d")
data_on.attrs["experiment_type"] = "nmr_spectrum"

data_on = dnp.remove_background(data_on)
data_on = dnp.apodize(data_on, lw=15)
data_on = dnp.fourier_transform(data_on)

# %%
# Plot Microwave On/Off DNP Spectra
# ---------------------------------
# First plot spectra individually

sampleTag = "10 mM TEMPO in Water"

dnp.plt.figure()
dnp.fancy_plot(data_on, title=sampleTag + ", MW On Spectrum")

dnp.plt.figure()
dnp.fancy_plot(data_off, title=sampleTag + ", MW Off Spectrum")
dnp.plt.show()

# %%
# Next plot both spectra in the same figure

dnp.plt.figure()
dnp.fancy_plot(data_on, xlim=[-20, 20])
dnp.fancy_plot(data_off*50, xlim=[-20, 20])
dnp.plt.title(sampleTag + ", MW ON/OFF(*50)")
dnp.plt.show()
