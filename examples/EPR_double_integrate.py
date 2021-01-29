"""
DNPLab EPR data tutorial
============================

This example demonstrates how to use DNPLab to load and double integrate EPR data.

"""
# Import DNPLab
import dnplab as dnp

# pro tip: shorten the syntax by aiming to specific functions if you would like. For example, get the 'load' function of the dnpImport module and use it without dnplab.dnpImport. in front.
from dnplab.dnpImport import load

# lets use some 1D topspin data. All data enters into the same object structure, so this example applies to any NMR format. EPR data are demonstrated in a separate example.
data = load("../data/bes3t/1D_CW.DTA", data_type="xenon")
# pro tip: if you do not shorten the syntax using the above tip, make sure to use:
#            data = dnplab.dnpImport.load(...
# pro tip: the 'data_type' argument is generally not needed, dnpImport autodetects the type.

# create a workspace, this is a place to store multiple objects if multiple processing methods are to be compared.
ws = dnp.create_workspace()
# add the data to the key 'raw'
ws.add("raw", data)
# copy the data to the processing buffer 'proc'
ws.copy("raw", "proc")
# raw data are now in both ws["proc"] and ws["raw"]. The following processing steps will modify ws["proc"] but leave ws["raw"] unotuched so that you can always return to the original data.

# pro tip: to see what attributes are available on import
print("BEFORE PROCESSING: " + str(ws["raw"].attrs.keys()))

# rename x axis to NMR style notation to more seamlessly interact with some functions
ws["proc"].rename("G", "f2")

# pro tip: see axis names with
print("raw axis names: " + str(ws["raw"].dims))
print("proc axis names: " + str(ws["proc"].dims))

# baseline correct with a zeroth order polynomial to remove the offset
dnp.dnpTools.baseline(ws, type="polynomial", order=0)

# pro tip: keep a copy of the processed spectrum before integrating
ws.copy("proc", "spec")

# double integrate
dnp.dnpTools.integrate(ws, type="double")

# access your processed spectrum as follows:
# field axis
x_axis = ws["raw"].coords["G"]
# spectrum
spectrum = ws["raw"].values
# first integral
first_integral = ws["proc"].attrs["first_integral"]
# double integral value
double_integral = ws["proc"].values

# print out the double integral value
print("Double integral = " + str(double_integral))

# pro tip: to see if any attributes have been added after processing
print("AFTER PROCESSING: " + str(ws["proc"].attrs.keys()))

# plot the EPR spectrum using dnpResults
dnp.dnpResults.plot(ws["raw"], label="EPR Spectrum")
dnp.dnpResults.legend()
dnp.dnpResults.show()

# pro tip: use matplotlib to plot arrays such as the attributes "first_integral" or "baseline" over the spectrum, in this case access the spectrum directly as '.values'
import matplotlib.pyplot as plt

plt.figure
plt.plot(ws["raw"].coords["G"], ws["spec"].values, label="corrected EPR Spectrum")
plt.plot(ws["raw"].coords["G"], ws["raw"].values, label="original EPR Spectrum")
plt.plot(ws["raw"].coords["G"], ws["proc"].attrs["baseline"], label="Baseline")
plt.plot(
    ws["raw"].coords["G"], ws["proc"].attrs["first_integral"], label="First Integral"
)
plt.legend()
plt.show()
