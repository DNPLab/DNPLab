import matplotlib.pyplot as _plt
import numpy as _np
from warnings import warn

from ..core.data import DNPData

dark_green = "#46812B"
light_green = "#67AE3E"
dark_grey = "#4D4D4F"
light_grey = "#A7A9AC"
orange = "#F37021"

_plt.rcParams["lines.linewidth"] = 1.5
_plt.rcParams["axes.prop_cycle"] = _plt.cycler(
    color=[orange, dark_green, light_green, dark_grey, light_grey]
)

show = _plt.show

# hand curated list of plotting arguments that are forwarded, note that this should probably be in a config file (refactoring needed)
forwarded_pyplot_plots = [
    "semilogy",
    "semilogx",
    "polar",
    "loglog",
    "scatter",
    "errorbar",
    "step",
]


def plot(data, *args, **kwargs):
    """Plot function for dnpdata object

    Args:
        data (DNPData): DNPData object for matplotlib plot function
        args: args for matplotlib plot function
        kwargs: kwargs for matplotlib plot function

        if any of semilogy, semilogx, polar, loglog, scatter, errorbar or step is in kwargs the argument will be evaluated with
        bool(). If this evaluates to True the corresponding matplotlib function is used instead of the standard plot

    Returns:
        Returns formated matplotlib plot.

    Example:

       Plotting a DNPData object:

       >>> dnp.plt.figure()
       >>> dnp.plot(data)
       >>> dnp.plt.show()

       # Plotting two DNPData objects (data1 and data2) on the same figure:

       >>> dnp.plt.figure()
       >>> dnp.plot(data1)
       >>> dnp.plot(data2)
       >>> dnp.plt.show()

       Plotting a DNPData object with some custom parameters:

       >>> dnp.plt.figure()
       >>> dnp.plot(data, 'k-', linewidth = 3.0, alpha = 0.5)
       >>> dnp.plt.show()

       Plotting a DNPData object with a semilogy plot (possible arguments: semilogy=1, semilogy=True, semilogy="True")
       Forwarded arguments: semilogy, semilogx, polar, loglog, scatter, errorbar or step
       The absolute value is taken to ensure that the y axis is always positive

       >>> dnp.plt.figure()
       >>> dnp.plot(np.abs(data), 'k-', linewidth = 3.0, alpha = 0.5, semilogy=1)
       >>> dnp.plt.show()

    """

    if "dim" in kwargs:
        dim = kwargs.pop("dim")
    else:
        dim = data.dims[0]

    coord = data.coords[dim]

    data.unfold(dim)

    # will try to plot various pyplot utility plot functions into same axis, the use should know what he does!
    # no unittest added, but only hand tested with semilogy and normal plot works as intended ni fancy_plot)
    use_default = True
    plot_function_list = []
    for k in forwarded_pyplot_plots:
        if bool(kwargs.pop(k, None)):
            plot_function_list.append(getattr(_plt, k))
            use_default = False
    for f in plot_function_list:
        f(coord, data.values.real, *args, **kwargs)
    if use_default:
        _plt.plot(coord, data.values.real, *args, **kwargs)
    data.fold()


def fancy_plot(data, xlim=[], title="", showPar=False, *args, **kwargs):
    """Streamline Plot function for dnpdata objects

    This function creates streamlined plots for NMR and EPR spectra. The type of the spectrum is detected from the attribute "experiment_type" of the dnpdata object. Currently the following types are implemented: nmr_spectrum, epr_spectrum, enhancements_P, and inversion_recovery. The function will automatically format the plot according to the "experiment_type" attribute.

    Args:
        data (DNPData): DNPData object with values to plot
        xlim (tuple): List of limit values for plotting function
        title (str): Plot title
        showPar (boolean): Toggle whether to show experiment parameters

    Returns:
        Returns formatted matplotlib plot.

    Example:

        Simply just plotting the dnpdata object:

        >>> dnp.fancy_plot(data)

        Plot EPR spectrum from 344 mT to 354 mT, show experimental parameters:

        >>> dnp.fancy_plot(data, xlim=[344, 354], title="EPR Spectrum", showPar=True)

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

    _plt.grid(True)
    _plt.title(title)

    if data.attrs["experiment_type"] == "nmr_spectrum":
        coord = data.coords[dim]
        data.unfold(dim)

        _plt.plot(coord, data.values.real, *args, **kwargs)
        _plt.xlabel("Chemical Shift $\delta$ (ppm)")
        _plt.ylabel("NMR Signal Intensity (a.u.)")

        _plt.xlim(max(coord), min(coord))

        if xlim != []:
            _plt.xlim(xlim[1], xlim[0])

        if showPar == True:
            parameterString = "Freq: " + str(round(data.attrs["nmr_frequency"], 4))

            box_style = dict(boxstyle="round", facecolor="white", alpha=0.25)
            xmin, xmax, ymin, ymax = _plt.axis()

            _plt.text(xmin * 0.95, ymax / 10, parameterString, bbox=box_style)

    elif data.attrs["experiment_type"] == "saturation_recovery":
        coord = data.coords[dim]
        data.unfold(dim)

        _plt.plot(
            coord, data.values.real, marker="o", fillstyle="none", *args, **kwargs
        )
        _plt.xlabel("Evolution Time T1 [s]")
        _plt.ylabel("Signal Intensity (a.u.)")

        if xlim != []:
            _plt.xlim(xlim[1], xlim[0])

        if title == "":
            _plt.title("Saturation Recovery")
        else:
            _plt.title(title)

    elif data.attrs["experiment_type"] == "polarization_buildup":
        coord = data.coords[dim]
        data.unfold(dim)

        _plt.plot(
            coord, data.values.real, marker="o", fillstyle="none", *args, **kwargs
        )
        _plt.xlabel("Contact Time t$_c$ [s]")
        _plt.ylabel("Signal Intensity (a.u.)")

        if xlim != []:
            _plt.xlim(xlim[1], xlim[0])

        if title == "":
            _plt.title("Polarization Build-Up")
        else:
            _plt.title(title)

    elif data.attrs["experiment_type"] == "dnp_enhancement_profile_f":
        coord = data.coords[dim]
        data.unfold(dim)

        _plt.plot(coord, data.values.real, *args, **kwargs)
        _plt.xlabel("Frequency [GHz]")
        _plt.ylabel("DNP Intensity (a.u.)")

        if xlim != []:
            _plt.xlim(xlim[1], xlim[0])

        if title == "":
            _plt.title("DNP Enhancement Profile")
        else:
            _plt.title(title)

    elif data.attrs["experiment_type"] == "eldor_profile":
        coord = data.coords[dim]
        data.unfold(dim)

        _plt.plot(coord, data.values.real, *args, **kwargs)
        _plt.xlabel("Frequency [GHz]")
        _plt.ylabel("ELDOR Intensity (a.u.)")

        if xlim != []:
            _plt.xlim(xlim[1], xlim[0])

        if title == "":
            _plt.title("ELDOR Spectrum")
        else:
            _plt.title(title)

    elif data.attrs["experiment_type"] == "echo_decay":
        coord = data.coords[dim]
        data.unfold(dim)

        _plt.plot(
            coord, data.values.real, marker="o", fillstyle="none", *args, **kwargs
        )
        _plt.xlabel("Decay Time [s]")
        _plt.ylabel("Signal Intensity (a.u.)")

        if xlim != []:
            _plt.xlim(xlim[1], xlim[0])

        if title == "":
            _plt.title("Echo Decay")
        else:
            _plt.title(title)

    elif data.attrs["experiment_type"] == "epr_spectrum":
        coord = data.coords[dim]
        data.unfold(dim)

        _plt.plot(coord, data.values.real, *args, **kwargs)
        _plt.xlabel("Magnetic Field $B_{0}$ (mT)")
        _plt.ylabel("EPR Signal Intensity (a.u.)")

        if xlim != []:
            _plt.xlim(xlim[0], xlim[1])

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
            xmin, xmax, ymin, ymax = _plt.axis()

            _plt.text(xmin * 1.001, ymin * 0.90, parameterString, bbox=box_style)

    elif (
        data.attrs["experiment_type"] == "enhancements_P"
        or data.attrs["experiment_type"] == "enhancements_PdBm"
    ):
        coord = data.coords[dim]
        data.unfold(dim)

        _plt.plot(
            coord, data.values.real, marker="o", fillstyle="none", *args, **kwargs
        )
        _plt.xlabel("Microwave Power (dBm)")
        _plt.ylabel("DNP Enhancement")

        if xlim != []:
            _plt.xlim(xlim[0], xlim[1])

        if title == "":
            _plt.title("DNP Enhancement Power Build-Up")
        else:
            _plt.title(title)

        # if showPar == True:

    elif data.attrs["experiment_type"] == "enhancements_PW":
        coord = data.coords[dim]
        data.unfold(dim)

        _plt.plot(
            coord,
            data.values.real,
            marker="o",
            fillstyle="none",
            *args,
            **kwargs
        )
        _plt.xlabel("Microwave Power (W)")
        _plt.ylabel("DNP Enhancement")

        if xlim != []:
            _plt.xlim(xlim[0], xlim[1])

        if title == "":
            _plt.title("DNP Enhancement Power Build-Up")
        else:
            _plt.title(title)

        # if showPar == True:

    elif data.attrs["experiment_type"] == "inversion_recovery":
        _plt.plot(
            data.coords["t1"],
            data.values.real,
            marker="o",
            fillstyle="none",
            *args,
            **kwargs
        )

        _plt.xlabel("Evolution Time T1 (s)")
        _plt.ylabel("NMR Amplitude [a.u.]")

        if xlim != []:
            _plt.xlim(xlim[0], xlim[1])

        if title == "":
            _plt.title("Inversion Recovery")
        else:
            _plt.title(title)

        # if showPar == True:

    else:
        plot(data, *args, **kwargs)

    data.fold()
