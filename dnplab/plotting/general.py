import matplotlib.pyplot as plt
import numpy as np

from ..core.data import DNPData

dark_green = "#46812B"
light_green = "#67AE3E"
dark_grey = "#4D4D4F"
light_grey = "#A7A9AC"
orange = "#F37021"

plt.rcParams["lines.linewidth"] = 1.5
plt.rcParams["axes.prop_cycle"] = plt.cycler(
    color=[orange, dark_green, light_green, dark_grey, light_grey]
)


def plot(data, *args, **kwargs):
    """Plot function for dnpdata object

    Args:
        data (dnpdata): dnpdata object for matplotlib plot function
        args: args for matplotlib plot function
        kwargs: kwargs for matplotlib plot function

    Example::

       # Plotting a dnpdata object
       dnp.dnpResults.plt.figure()
       dnp.dnpResults.plot(data)
       dnp.dnpResults.plt.show()

       # Plotting a workspace (dnpdata_collection)
       dnp.dnpResults.plt.figure()
       dnp.dnpResults.plot(ws['proc'])
       dnp.dnpResults.plt.show()

       # Plotting two curves on the same figure
       dnp.dnpResults.plt.figure()
       dnp.dnpResults.plot(ws['proc1'])
       dnp.dnpResults.plot(ws['proc2'])
       dnp.dnpResults.plt.show()

       # Plotting with some custom parameters
       dnp.dnpResults.plt.figure()
       dnp.dnpResults.plot(ws['proc'], 'k-', linewidth = 3.0, alpha = 0.5)
       dnp.dnpResults.plt.show()

    """

    if "dim" in kwargs:
        dim = kwargs.pop("dim")
    else:
        dim = data.dims[0]

    coord = data.coords[dim]

    data.unfold(dim)

    plt.plot(coord, data.values.real, *args, **kwargs)
    data.fold()


def dnplabplot(data, xlim = [], title = '', *args, **kwargs):
    """Streamline Plot function for dnpdata objects

    This function creates streamlined plots for NMR and EPR spectra. The type of the spectrum is automatically detected and the axis are formated accordingly.

    Args:

    Returns:

    Example:
    """

    if "dim" in kwargs:
        dim = kwargs.pop("dim")
    else:
        dim = data.dims[0]

    coord = data.coords[dim]
    data.unfold(dim)

    if dim == "f2":
        plt.plot(coord, data.values.real, *args, **kwargs)
        plt.grid(True)

        plt.xlabel("Chemical Shift $\delta$ (ppm)")
        plt.ylabel("NMR Signal Intensity (a.u.)")
        plt.title(title)

        plt.xlim(max(coord), min(coord))

        if xlim != []:

            plt.xlim(xlim[1],xlim[0])

    elif dim == "B0":
        plt.plot(coord, data.values.real, *args, **kwargs)
        plt.grid(True)
        plt.xlabel("Magnetic Field $B_{0}$ (mT)")
        plt.ylabel("EPR Signal Intensity (a.u.)")
        plt.title(title)

        if xlim != []:

            plt.xlim(xlim[0],xlim[1])

    else:
        print("Spectrum type unknown")

    data.fold()