import matplotlib.pyplot as plt
import numpy as np
from warnings import warn

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

show = plt.show


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


def fancy_plot(data, xlim=[], title="", showPar=False, *args, **kwargs):
    """Streamline Plot function for dnpdata objects

    This function creates streamlined plots for NMR and EPR spectra. The type of the spectrum is detected from the attribute "experiment_type" of the dnpdata object. Currently the following types are implemented: nmr_spectrum, epr_spectrum, enhancements_P, and inversion_recovery.

    Args:
        data (dnpdata): dnpdata object with values to plot
        xlim: list of limit values for plotting function
        title: string containing plot title
        showPar: boolean, toggle whether to show experiment parameters

    Returns:
        Returns formated matplotlib plot.

    Example::

        # Simply just plotting the dnpdata object
        dnp.fancy_plot(data)

        # Plot EPR spectrum from 344 mT to 354 mT, show experimental parameters
        dnp.fancy_plot(data, xlim=[344, 354], title="EPR Spectrum", showPar=True)

    """

    if "experiment_type" not in data.attrs:
        warn("experiment_type not defined in data.attrs, falling back to plot function")
        plot(data, *args, **kwargs)
        return
    elif data.attrs["experiment_type"] == None:
        warn("experiment_type is None, falling back to plot function")
        plot(data, *args, **kwargs)
        return

    if "dim" in kwargs:
        dim = kwargs.pop("dim")
    else:
        dim = data.dims[0]

    plt.grid(True)
    plt.title(title)

    if data.attrs["experiment_type"] == "nmr_spectrum":
        coord = data.coords[dim]
        data.unfold(dim)

        plt.plot(coord, data.values.real, *args, **kwargs)
        plt.xlabel("Chemical Shift $\delta$ (ppm)")
        plt.ylabel("NMR Signal Intensity (a.u.)")

        plt.xlim(max(coord), min(coord))

        if xlim != []:

            plt.xlim(xlim[1], xlim[0])

        if showPar == True:

            parameterString = "Freq: " + str(round(data.attrs["nmr_frequency"], 4))

            box_style = dict(boxstyle="round", facecolor="white", alpha=0.25)
            xmin, xmax, ymin, ymax = plt.axis()

            plt.text(xmin * 0.95, ymax / 10, parameterString, bbox=box_style)

    elif data.attrs["experiment_type"] == "epr_spectrum":
        coord = data.coords[dim]
        data.unfold(dim)

        plt.plot(coord, data.values.real, *args, **kwargs)
        plt.xlabel("Magnetic Field $B_{0}$ (mT)")
        plt.ylabel("EPR Signal Intensity (a.u.)")

        if xlim != []:
            plt.xlim(xlim[0], xlim[1])

        if showPar == True:
            SW = coord[-1] - coord[0]

            parameterString = (
                "MF: "
                + str(round(data.attrs["frequency"], 4))
                + "\nMP: "
                + str(round(data.attrs["power"], 3))
                + "\nCF: "
                + str(round(data.attrs["center_field"] / 10, 2))
                + "\nSW: "
                + str(round(SW, 2))
                + "\nMA: "
                + str(round(data.attrs["modulation_amplitude"], 2))
                + "\nNS: "
                + str(data.attrs["nscans"])
                + "\nTM: "
                + str(round(data.attrs["temperature"], 1))
            )

            box_style = dict(boxstyle="round", facecolor="white", alpha=0.25)
            xmin, xmax, ymin, ymax = plt.axis()

            plt.text(xmin * 1.001, ymin * 0.90, parameterString, bbox=box_style)

    elif data.attrs["experiment_type"] == "enhancements_P":
        coord = data.coords[dim]
        data.unfold(dim)

        plt.plot(coord, data.values.real, marker="o", fillstyle="none", *args, **kwargs)
        plt.xlabel("Microwave Power (dBm)")
        plt.ylabel("DNP Enhancement")

        if xlim != []:
            plt.xlim(xlim[0], xlim[1])

        # if showPar == True:

    elif data.attrs["experiment_type"] == "inversion_recovery":
        plt.plot(
            data.coords["t1"],
            data.values.real,
            marker="o",
            fillstyle="none",
            *args,
            **kwargs
        )

        plt.xlabel("Evolution Time T1 (s)")
        plt.ylabel("NMR Amplitude [a.u.]")

        if xlim != []:
            plt.xlim(xlim[0], xlim[1])

        # if showPar == True:

    else:

        plot(data)

    data.fold()
