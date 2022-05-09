# %%
"""
.. _plot_05_apodization:

===========
Apodization
===========

This example how to use the apodization (windowing functions) in DNPLab.

"""
# %%
# First we have to prepare the Python environment by importing DNPLab and load an example spectrum.

import dnplab as dnp
import matplotlib.pylab as plt

import numpy as np

# data = dnp.load("../../data/topspin/1")
data_raw = dnp.load("../../data/prospa/toluene_10mM_Tempone/41/data.1d")
data_raw.attrs["experiment_type"] = "nmr_spectrum"

# %%
# The spectral data is now stored in the data object. DNPLab uses object oriented programming, therefore, data is not just a variable. It is a dnpData object. Next, we will perform some basic processing, methods that are pretty standard in NMR spectroscopy. For this just pass the dnpData object from one function to the next.

# data = dnp.remove_background(data)
data = dnp.apodize(data_raw, lw=1)
data = dnp.fourier_transform(data)

# %%
# Reference spectrum. This will be more convenient in future versions of DNPLab

ref = 3.55
coord = data.coords["f2"]
data.coords["f2"] = coord + ref

# data = dnp.autophase(data, method="manual", phase=-1 * dnp.constants.pi / 2)

# %%
# Or by using the autophase function to automatically search for and apply the best phase angle

data = dnp.autophase(data, method="manual", phase=175 * dnp.constants.pi / 180)
# dnp.fancy_plot(data, xlim = [-15, 25])
# dnp.plt.show()


# %%
# To change the linewidth to 10 Hz do this

data_exp = dnp.apodize(data_raw, kind="exponential",lw=10)
data_exp = dnp.fourier_transform(data_exp)

ref = 3.55
coord = data_exp.coords["f2"]
data_exp.coords["f2"] = coord + ref
data_exp = dnp.autophase(data_exp, method="manual", phase=175 * dnp.constants.pi / 180)


data_gauss = dnp.apodize(data_raw, kind="gaussian", lw=10)
data_gauss = dnp.fourier_transform(data_gauss)

ref = 3.55
coord = data_gauss.coords["f2"]
data_gauss.coords["f2"] = coord + ref
data_gauss = dnp.autophase(
    data_gauss, method="manual", phase=175 * dnp.constants.pi / 180
)

data_gauss.values = data_gauss.values + 1000

# dnp.fancy_plot(data_exp, xlim=[-15, 25], label="Exponential")
# dnp.fancy_plot(data_gauss, xlim=[-15, 25], label="Gaussian")
# dnp.plt.legend()

# dnp.plt.show()




x = np.linspace(0,100,100)
lw = 10

y_exp = dnp.window.exponential(x,lw)
y_gauss = dnp.window.gaussian(x,lw)

dnp.plt.plot(x,y_exp,label = "Exponential")
dnp.plt.plot(x,y_gauss,label="Gaussian")
dnp.plt.legend()
dnp.plt.show()



