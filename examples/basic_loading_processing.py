"""
DNPLab basics tutorial
============================

This example demonstrates how to import and process data. The user also learns how to access the various components of the data object and generally how to interact with the functions.

"""
# Import standard modules if needed
import numpy as np

# Import DNPLab
import dnplab as dnp

# pro tip: shorten the syntax by aiming to specific functions if you would like. For example,
# get the 'load' function of the dnpImport module and use it without dnplab.dnpImport.
from dnplab.dnpImport import load

# lets use some 1D topspin data. All data enters into the same object structure, so this
# example applies to any NMR format. EPR data are demonstrated in a separate example.
data = load("../data/topspin/5", data_type="topspin")
# pro tip: if you do not shorten the syntax using the above tip, make sure to use:
#            data = dnplab.dnpImport.load(...
# pro tip: the 'data_type' argument is generally not needed, dnpImport autodetects the type

# create a workspace, this is a place to store multiple objects if multiple processing
# methods are to be compared.
ws = dnp.create_workspace()
# add the data to the key 'raw'
ws.add("raw", data)
# copy the data to the processing buffer 'proc'
ws.copy("raw", "proc")
# raw data are now in both ws["proc"] and ws["raw"]. The following processing steps will
# modify ws["proc"] but leave ws["raw"] unotuched so that you can always return to the
# original data.

# lets perform some basic processing. simply pass the entire workspace, 'ws', to the
# functions and the ws["proc"] will be manipulated. For example:
dnp.dnpNMR.remove_offset(ws)
dnp.dnpNMR.window(ws, linewidth=10)
dnp.dnpNMR.fourier_transform(ws, zero_fill_factor=2)
# this procedure removes the DC offset, multiplies by a window function, and fourier
# transforms.

# you can manually phase correct using the 'autophase' function using
dnp.dnpNMR.autophase(ws, method="manual", phase=np.pi / 2)
# or use autophase to search for and apply the best phase angle
dnp.dnpNMR.autophase(ws)
# pro tip: default order is zeroth, first order is also available and is demonstrated
# below.

# access your processed spectrum as follows:
# ppm axis
x_axis = ws["proc"].coords["f2"]
# spectrum
spectrum = ws["proc"].values
# dimensions by default are: before FT, direct is "t2"; indirect is "t1";
#                            after FT, direct is "f2"; indirect is still "t1";
# rename dimensions at anytime with, for example,
ws["proc"].rename("f2", "ppm")

# access data attributes such as NMR frequency as follows:
nmr_frequency = ws["proc"].attrs["nmr_frequency"]
# pro tip: to see what attributes are available use
print(ws["proc"].attrs.keys())

# plot the frequency spectrum, imaginary and real
dnp.dnpResults.plot(ws["proc"].imag, label="Imaginary")
dnp.dnpResults.plot(ws["proc"].real, label="Real")
dnp.dnpResults.legend()
dnp.dnpResults.show()

# now suppose you have 2D inversion recovery data, import using
data = load("../data/topspin/304")
ws_2D = dnp.create_workspace()
ws_2D.add("raw", data)
ws_2D.copy("raw", "proc")

# process as above
dnp.dnpNMR.remove_offset(ws_2D)
dnp.dnpNMR.window(ws_2D, linewidth=10)
dnp.dnpNMR.fourier_transform(ws_2D, zero_fill_factor=2)
# this time use a first order phase correction where the phase angle changes ('delta')
# linearly by pi/16 over the spectral range and pivots 100 points into the spectrum.
dnp.dnpNMR.autophase(ws_2D, method="arctan", order="first", pivot=100, delta=np.pi / 16)

# now align the spectra
dnp.dnpNMR.align(ws_2D)
# store a copy of the data before they are reduced to integrals
ws_2D.copy("proc", "spec")
# now integrate
dnp.dnpTools.integrate(ws_2D)
# pro tip: by default the entire spectrum is integrated, be sure to optimize the width
# of the center and width of the integration window. For example,

# so lets integrate only over a certain range instead. First return to the data as
# it was before integrating by copying the saved 'spec'' back to 'proc'
ws_2D.copy("spec", "proc")
# integrate over a 100ppm range centered at 10ppm
dnp.dnpTools.integrate(ws_2D, integrate_center=10, integrate_width=100)

# import the exponential fitting function with shortened syntax
from dnplab.dnpFit import exponential_fit as fit

# fit the integrals to an inversion recovery function to extract T1 using fit type="T1".
# Fit data are added to the object as the key 'fit'.
fit(ws_2D, type="T1")

# access the fit data as follows
# indirect axis
x_axis = ws_2D["fit"].coords["t1"]
# fit
fit_curve = ws_2D["fit"].values
# T1 value
T1 = ws_2D["fit"].attrs["T1"]
# standard deviation in T1 calculated using the jacobian->covariance->stdd method
T1_stdd = ws_2D["fit"].attrs["T1_stdd"]

# plot the integrals and T1 fit
dnp.dnpResults.figure()
dnp.dnpResults.plot(ws_2D["proc"].imag, marker="o", linestyle="none", label="Imaginary")
dnp.dnpResults.plot(ws_2D["proc"].real, marker="o", linestyle="none", label="Real")
dnp.dnpResults.plot(ws_2D["fit"], label="Fit")
dnp.dnpResults.legend()
dnp.dnpResults.show()
