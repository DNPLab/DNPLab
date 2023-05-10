import matplotlib.pyplot as _plt
import numpy as _np
from warnings import warn

from ..core.data import DNPData
from ..__init__ import DNPLAB_CONFIG
from ..__init__ import _get_dnp_config

"""
Standard dnplab colors

dark_green = DNPLAB_CONFIG.get('COLORS','dark_green')
light_green = DNPLAB_CONFIG.get('COLORS','light_green')
dark_grey = DNPLAB_CONFIG.get('COLORS','dark_grey')
light_grey = DNPLAB_CONFIG.get('COLORS','light_grey')
orange=DNPLAB_CONFIG.get('COLORS','orange')
"""

FANCYPLOT_CONFIG=_get_dnp_config(DNPLAB_CONFIG.get('CONFIGNAMES','FANCYPLOT_CONFIG'))
# hand curated list of plotting arguments that are forwarded, from config file
forwarded_pyplot_plots=DNPLAB_CONFIG.getlist('PLOTTING','forwarded_pyplot_plots')


cycler_list = [ DNPLAB_CONFIG.get('COLORS',color_key) for color_key in DNPLAB_CONFIG['COLORS'].keys() ]
_plt.rcParams["lines.linewidth"] = 1.5
_plt.rcParams["axes.prop_cycle"] = _plt.cycler(
    color= cycler_list
)

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

    _plt.grid(True)
    _plt.title(title)

    if data.attrs["experiment_type"] == "nmr_spectrum":
        if "dim" in kwargs:
            dim = kwargs.pop("dim")
        else:
            dim = data.dims[0]
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

    if data.attrs["experiment_type"] in FANCYPLOT_CONFIG.sections():
        exp_type=data.attrs["experiment_type"]
        get_key = lambda x,fallback=None: FANCYPLOT_CONFIG.get(exp_type,x,fallback=fallback)
        get_float_key = lambda x,fallback=1: FANCYPLOT_CONFIG.getfloat(exp_type,x,fallback=fallback)

        dim = kwargs.pop("dim",FANCYPLOT_CONFIG.get(exp_type,"dim",fallback=data.dims[0]))
        coord = data.coords[dim] * get_float_key('coord_scaling')
        data.unfold(dim)

        plt_config_kwargs = {key.lstrip('__'):val for key,val in FANCYPLOT_CONFIG[exp_type].items() if key.startswith('__') }
        plt_config_kwargs.update(kwargs) #calling values take precedence over config values

        _plt.plot(
            coord, data.values.real * get_float_key('value_scaling') , *args, **plt_config_kwargs
        )
        _plt.xlabel(get_key('xlabel'))
        _plt.ylabel(get_key("ylabel"))

        if xlim != []:
            _plt.xlim(xlim[1], xlim[0])

        if title == "":
            _plt.title(get_key("title"))
        else:
            _plt.title(title)

        if  FANCYPLOT_CONFIG.getboolean(exp_type,'showPar',fallback=False) or showPar:
            SW = coord[-1] - coord[0]
            # alternative: if "showPar" in FANCYPLOT_CONFIG[exp_type].keys():
            # (MF,frequency,4);(MP,power,3);....
            # attrs=[k.strip("(").strip(")").split(",") for k in FANCYPLOT_CONFIG[exp_type]["showPar"].split(';')]
            # prmString=''
            # for k in attrs:
            #   try:
            #      label,attribute, round_value = k
            #   except ValueError as ve:
            #       warn("could not unpack attribute tupel {0} into label, attribute and round_value, skipping this entry!".format(k))
            #       continue
            #   try:
            #       try:
            #          prmString+=label.strip().strip(":")+": "+round(data.attrs[attribute],int(round_value))+"\n"
            #       except KeyError as e:
            #           warn("Attribute {0} not in data.attributes, skipping this entry!")
            #   except ValueError as ve:
            #       warn("Could not convert {0} to integer, skipping this entry because of error {1}".format(round_value,ve))

            try:
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
            except KeyError as e:
                warn('Trying to show parameters but error {0} occured. Are you sure that the attributes are in the data.attrs?\nParameterstring is empty!'.format(e))
                parameterString=''

            box_style = dict(boxstyle="round", facecolor="white", alpha=0.25)
            xmin, xmax, ymin, ymax = _plt.axis()

            _plt.text(xmin * 1.001, ymin * 0.90, parameterString, bbox=box_style)
        ax=_plt.gca()
        fig=_plt.gcf()
        for key in FANCYPLOT_CONFIG[exp_type].keys():
            value=FANCYPLOT_CONFIG[exp_type][key]
            if key.startswith('ax.'):
                ax_key=key.lstrip('ax.')
                try:
                    getattr(ax,ax_key)(value)
                except ValueError as e:
                    warn("Could not set ax attribute {0} to string value {1}, skipping this option! (ValueError: {2})".format(ax_key,value,e))
            if key.startswith('fig.'):
                fig_key=key.lstrip('fig.')
                try:
                    getattr(fig,fig_key)(value)
                except ValueError as e:
                    warn("Could not set figure attribute {0} to string value {1}, skipping this option! (ValueError: {2})".format(fig_key,value,e))

    else:
        plot(data, *args, **kwargs)

    data.fold()
