# %%
"""
.. _02_load_1d_spectrum_kea:

===================================
02 - Load 1D NMR spectra (Kea Data)
===================================

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
data_off = dnp.load("../data/prospa/10mM_TEMPO_Water/1Pulse_20200929/35/data.1d")

data_off = dnp.remove_background(data_off)
data_off = dnp.apodize(data_off, lw = 15)
data_off = dnp.fourier_transform(data_off)

# %%
# Load and Process ON-Signal
# ----------------------------
# Importing the on-signal involves the same steps as importing the off-signal. Once processed the data is copied to the results buffer 'onSignal'.

########## ON Signal (P = 2 W) ##########
data_on = dnp.load("../data/prospa/10mM_TEMPO_Water/1Pulse_20200929/51/data.1d")

data_on = dnp.remove_background(data_on)
data_on = dnp.apodize(data_on, lw = 15)
data_on = dnp.fourier_transform(data_on)

# %%
# Plot Microwave On/Off DNP Spectra
# ---------------------------------
# First plot individual spectra

sampleTag = "10 mM TEMPO in Water, "

# dnp.dnplabplot(data_on, title = sampleTag + "MW ON Spectrum")
# dnp.dnplabplot(data_off, title = sampleTag + "MW ON Spectrum")

# %%
# Next plot both spectra in the same figure

dnp.plt.figure()
dnp.plt.plot(data_on.coords["f2"], data_on.values.real, label = "MW On Signal (2 W)")
dnp.plt.plot(data_on.coords["f2"], data_off.values.real * 10 - 100, label = "MW Off Signal x 10")
dnp.plt.xlim([30, -30])
dnp.plt.xlabel("Chemical Shift [ppm]")
dnp.plt.ylabel("Signal Amplitude [a.u.]")
dnp.plt.title("DNP On/Off Signal, 10 mM TEMPO in Water")
dnp.plt.legend()
dnp.plt.grid(True)
dnp.plt.show()