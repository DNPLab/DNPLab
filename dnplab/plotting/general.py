import matplotlib.pyplot as _plt
import numpy as _np
from warnings import warn

from ..core.data import DNPData
from ..config.config import DNPLAB_CONFIG

"""
Standard dnplab colors

dark_green = DNPLAB_CONFIG.get('COLORS','dark_green')
light_green = DNPLAB_CONFIG.get('COLORS','light_green')
dark_grey = DNPLAB_CONFIG.get('COLORS','dark_grey')
light_grey = DNPLAB_CONFIG.get('COLORS','light_grey')
orange=DNPLAB_CONFIG.get('COLORS','orange')
"""

# hand curated list of plotting arguments that are forwarded, from config file
_forwarded_pyplot_plots = DNPLAB_CONFIG.getlist("PLOTTING", "forwarded_pyplot_plots")


_cycler_list = [
    DNPLAB_CONFIG.get("COLORS", color_key)
    for color_key in DNPLAB_CONFIG["COLORS"].keys()
]
_plt.rcParams["lines.linewidth"] = 1.5
_plt.rcParams["axes.prop_cycle"] = _plt.cycler(color=_cycler_list)

# As discussed: for now change the rcParams, we do not use a temporary context - the values are stored in the dnplab config
_plt.rcParams["font.family"] = "Arial"
_plt.rcParams["pdf.fonttype"] = 42
_plt.rcParams["font.size"] = 14

# default context
# _dnpContext = {"font.family": "Arial", "pdf.fonttype": 42, "font.size":14}

show = _plt.show


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

    # default context

    # will try to plot various pyplot utility plot functions into same axis, the use should know what he does!
    # no unittest added, but only hand tested with semilogy and normal plot works as intended ni fancy_plot)
    use_default = True
    plot_function_list = []
    for k in _forwarded_pyplot_plots:
        k = k.strip()
        if bool(kwargs.pop(k, None)):
            plot_function_list.append(getattr(_plt, k))
            use_default = False
            plot_return = []
    for f in plot_function_list:
        plot_return.append(f(coord, data.values.real, *args, **kwargs))
    if use_default:
        plot_return = _plt.plot(coord, data.values.real, *args, **kwargs)
    fontsize_xlabel = DNPLAB_CONFIG.getint("PLOTTING", "fontsize_xlabel")
    _plt.xlabel(dim, fontsize=fontsize_xlabel)
    _plt.ylabel("", fontsize=fontsize_xlabel)
    data.fold()

    return plot_return


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

    _plt.grid(True)
    _plt.title(title)

    fancyplot_possiblesections = list(DNPLAB_CONFIG.sections())
    fancyplot_label = DNPLAB_CONFIG.get(
        "PLOTTING", "fancyplot_label", fallback="FANCY_PLOT"
    )
    fancyplot_sections = [
        k.strip(fancyplot_label).strip(":")
        for k in fancyplot_possiblesections
        if k.startswith(fancyplot_label)
    ]

    if data.attrs["experiment_type"] == "nmr_spectrum":
        if "dim" in kwargs:
            dim = kwargs.pop("dim")
        else:
            dim = data.dims[0]
        coord = data.coords[dim]
        data.unfold(dim)

        plot_return = _plt.plot(coord, data.values.real, *args, **kwargs)
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

        data.fold()

    elif data.attrs["experiment_type"] in fancyplot_sections:
        exp_type = fancyplot_label + ":" + data.attrs["experiment_type"]
        get_key = lambda x, fallback=None: DNPLAB_CONFIG.get(
            exp_type, x, fallback=fallback
        )
        get_float_key = lambda x, fallback=1: DNPLAB_CONFIG.getfloat(
            exp_type, x, fallback=fallback
        )

        dim = kwargs.pop(
            "dim", DNPLAB_CONFIG.get(exp_type, "dim", fallback=data.dims[0])
        )
        coord = data.coords[dim] * get_float_key("coord_scaling")
        data.unfold(dim)
        plt_config_kwargs = {
            key.lstrip("__"): val
            for key, val in DNPLAB_CONFIG[exp_type].items()
            if key.startswith("__")
        }
        plt_config_kwargs.update(
            kwargs
        )  # calling values take precedence over config values

        plot_return = _plt.plot(
            coord,
            data.values.real * get_float_key("value_scaling"),
            *args,
            **plt_config_kwargs
        )

        if xlim != []:
            _plt.xlim(xlim[0], xlim[1])

        ax = _plt.gca()
        fig = _plt.gcf()
        for key in DNPLAB_CONFIG[exp_type].keys():
            if key.startswith("ax.") or key.startswith("fig."):
                args, kwargs = DNPLAB_CONFIG.getargs_kwargs(exp_type, key)
                prm_key = key.lstrip("ax.").lstrip("fig.")
                try:
                    if key.startswith("ax."):
                        getattr(ax, prm_key)(*args, **kwargs)
                    else:
                        getattr(fig, prm_key)(*args, **kwargs)
                except ValueError as e:
                    warn(
                        "Could not set ax/fig attribute {0} to string value {1}, skipping this option! (ValueError: {2})".format(
                            prm_key, (args, kwargs), e
                        )
                    )

        if title != "":
            _plt.title(title)

        if showPar:
            prmString = ""
            keylist = list(DNPLAB_CONFIG[exp_type].keys())
            attrs_tpl = [
                (k.lstrip("showpar_"), k)
                for k in keylist
                if (k.startswith("showpar_") and (not k.endswith("_scaling")))
            ]
            for attr, key in attrs_tpl:
                try:
                    scaling = DNPLAB_CONFIG.getfloat(
                        exp_type, key + "_scaling", fallback=1
                    )
                    prmString += DNPLAB_CONFIG[exp_type][key].format(
                        data.attrs[attr] * scaling
                    )
                    prmString.strip()
                    if prmString[-1] != "\n":
                        prmString += "\n"
                except KeyError:
                    warn(
                        "Attribute {0} not in data.attributes, skipping this entry!".format(
                            attr
                        )
                    )

            SW = coord[-1] - coord[0]
            prmString += "SW: " + str(round(SW, 2))
            box_style = dict(boxstyle="round", facecolor="white", alpha=0.25)
            xmin, xmax, ymin, ymax = _plt.axis()

            _plt.text(xmin * 1.001, ymin * 0.90, prmString, bbox=box_style)

        data.fold()

    else:
        plot_return = plot(data, *args, **kwargs)

    _plt.tight_layout()

    return plot_return
