import matplotlib.pyplot as plt
import numpy as np

from ..core.data import DNPData

dark_green = "#46812B"
light_green = "#67AE3E"
dark_grey = "#4D4D4F"
light_grey = "#A7A9AC"
orange = "#F37021"

plt.rcParams["figure.figsize"] = [7, 5]
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["pdf.fonttype"] = 42

plt.rcParams["lines.linewidth"] = 1.5
plt.rcParams["font.size"] = 18
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


        import matplotlib.pylab as plt

        # Plotting a dnpdata object
        plt.figure()
        dnp.plot(data)
        plt.show()

        # Plotting with some custom parameters
        plt.figure()
        dnp.plot(data, 'k-', linewidth = 3.0, alpha = 0.5)
        plt.show()

    """

    if "dim" in kwargs:
        dim = kwargs.pop("dim")
    else:
        dim = data.dims[0]

    coord = data.coords[dim]

    data.unfold(dim)

    plt.plot(coord, data.values.real, *args, **kwargs)
    plt.xlabel(dim)
    data.fold()
