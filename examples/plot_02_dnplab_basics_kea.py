# %%
"""
.. _02_dnplab_baseics_kea_data:

=============================
02 - DNPLab Basics (Kea Data)
=============================

In this example we demonstrate how to import two NMR spectra, one recorded with a microwave power of 0 W (off-signal) and one with a microwave power of 2 W (on-signal). The spectra are recorded using a Magritek Kea system.

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
#
# The next section demonstrates how a FID is imported into DNPLab (lines 35) and how a workspace is created (line 37). Many functions within DNPLab work by default on the processing buffer ('proc') and leave the raw data untouched. After importing, the data is first copied to the 'rawOff' buffer (line 38) and then the 'rawOff' buffer is copied to the 'proc' buffer (line 39).
# Afterwards the FID is processed, first by removing any DC offset (line 41), then a 15 Hz linewidth is applied (line 42) before performing the Fourier transformation (line 43). At the end the processed data in the 'proc' buffer is copied to a results buffer called 'offSignal' 9line 45).

########## OFF Signal (P = 0 W) ##########
data = dnp.dnpImport.load("../data/prospa/10mM_TEMPO_Water/1Pulse_20200929/35/data.1d")

ws = dnp.create_workspace()
ws.add("rawOff", data)
ws.copy("rawOff", "proc")

dnp.dnpNMR.remove_offset(ws)
dnp.dnpNMR.window(ws, linewidth=15)
dnp.dnpNMR.fourier_transform(ws, zero_fill_factor=2)

ws.copy("proc", "offSignal")


# %%
# Load and Process ON-Signal
# ----------------------------
#
# Importing the on-signal involves the same steps as importing the off-signal. Once processed the data is copied to the results buffer 'onSignal'.

########## ON Signal (P = 2 W) ##########
data = dnp.dnpImport.load("../data/prospa/10mM_TEMPO_Water/1Pulse_20200929/51/data.1d")

ws.add("rawOn", data)
ws.copy("rawOn")

dnp.dnpNMR.remove_offset(ws)
dnp.dnpNMR.window(ws, linewidth=15)
dnp.dnpNMR.fourier_transform(ws, zero_fill_factor=2)

ws.copy("proc", "onSignal")


# %%
# Create a Figure and Plot On/Off Spectra
# ---------------------------------------
#
#

########## Plot Spectra ##########
dnp.dnpResults.figure()
dnp.dnpResults.plot(ws["offSignal"].real * 10 - 100, label="Off Signal x 10")
dnp.dnpResults.plot(ws["onSignal"].real, label="On Signal")
dnp.dnpResults.xlim([30, -30])
dnp.dnpResults.plt.xlabel("Chemical Shift [ppm]")
dnp.dnpResults.plt.ylabel("Signal Amplitude [a.u.]")
dnp.dnpResults.plt.title("DNP On/Off Signal, 10 mM TEMPO in Water")
dnp.dnpResults.plt.legend()
dnp.dnpResults.plt.grid(True)
dnp.dnpResults.show()
