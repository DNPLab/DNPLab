import matplotlib.pyplot as plt
import numpy as _np


def stack(data, *args, offset=None, **kwargs):
    """Stack Plot for 2D data

    Args:
        data (dnpdata): dnpdata object for matplotlib plot function
        args: args for matplotlib plot function
        offset: Value to offset each spectra, by default maximum of absolute value
        kwargs: kwargs for matplotlib plot function

    Examples:

       dnp.stack(data)

    """
    waterfall(data, 0, offset, *args, **kwargs)


def waterfall(data, dx, dy, *args, **kwargs):
    """Waterfall plot for 2d data

    Args:
        data (dnpData): 2d Data object for waterfall plot
        dx (float, int): x-increment for each line
        dy (float, int): y-increment for each line

    Examples:

       dnp.waterfall(data)

    """

    coord = data.coords[0]
    dim = data.dims[0]

    for ix in range(data.coords[data.dims[1]].size):
        plt.plot(
            coord + (ix * dx),
            data[data.dims[1], ix].values.ravel() + (ix * dy),
            *args,
            **kwargs,
            zorder=-1 * ix
        )
        plt.fill_between(
            coord + (ix * dx),
            data[data.dims[1], ix].values.ravel() + (ix * dy),
            ix * dy,
            facecolor="w",
            edgecolor="None",
            zorder=-1 * ix,
        )
