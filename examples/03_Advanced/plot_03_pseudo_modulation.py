# %%
"""
.. _plot_03_pseudo_modulation:

================================
Pseudo Modulation of EPR Spectra
================================

Typically, cw EPR spectra are recorded and shown as its first derivative, because the spectrum is detected using a lock-in amplifier. In contrast, a echo-detected field sweep spectrum is recorded as the absorption spectrum. To simulated the effect of the lock-in detection on an absorption spectrum to compare spectra, the data can be pseudo-modulated to calculate the derivative of the spectrum and to filter out noise. The same attention needs to be paid to the modulation amplitude as in an actual cw experiment. The spectrum can easily be overmodulated if the value of the modulation amplitude is too high.

The procedure was developed in the lab of Jim Hyde and is described here:

* Hyde, J., M. Pasenkiewicz-Gierula, A. Jesmanowicz, and W. Antholine. “Pseudo Field Modulation in EPR Spectroscopy.” Applied Magnetic Resonance 1 (1990): 483–96.

"""
# %%
# Load 1D EPR Spectrum
# --------------------
# For this example we will start by loading a echo-detected field sweep spectrum of a nitroxide radical. The spectrum was recorded using SpecMan4EPR using a constant magnetic field sweep. The field is recorded at each point and fitted to a 3rd order polynome to generate the field axis. Many of these steps are specific to the data set and may not required in your case.

import dnplab as dnp
import numpy as np

file_name_path = "../../data/specman/Nitroxide_Q_Band.d01"
data = dnp.load(file_name_path)

data = dnp.left_shift(data,'x0',1)

field = data['x1',2].values.squeeze()
p = np.polyfit(data.coords['x0'],field, deg = 3)
field_fit = np.polyval(p, data.coords['x0'])

data = data['x1',0].sum('x1') + 1j * data['x1',1].sum('x1')
data = dnp.update_axis(data, start_stop = [field_fit[0], field_fit[-1]], dim = 0, new_dims = "B0")

# %%
# The imported EPR spectrum is shown in the figure below. It is an echo-detected field-sweep EPR spectrum of a nitroxide bi-radical, recorded at Q-Band frequencies.

dnp.plt.figure()
dnp.fancy_plot(data)
dnp.plt.show()

# %%
# To calculated the pseudo-modulated spectrum of the absorption spectrum use the DNPLab function ``pseudo_modulation``. 

data_mod = dnp.pseudo_modulation(data, modulation_amplitude = 0.001)

# %%
# The field axis of this example spectrum is given in (T). Therefore, the unit for the modulation_amplitude also needs to be given in (T). In the above example, the pseudo-modulation spectrum is calculated using a modulation amplitude of 1 mT (10 G).

dnp.plt.figure()
dnp.fancy_plot(data_mod)
dnp.plt.show()

# %%
# Just like in a real cw experiment, the spectrum will show too broad lines, if the modulation amplitude is too large. The spectrum will be overmodulated. Below, the effect is shown for a modulation amplitude of 5 mT (50 G).

data_mod = dnp.pseudo_modulation(data, modulation_amplitude = 0.005)

dnp.plt.figure()
dnp.fancy_plot(data_mod)
dnp.plt.show()