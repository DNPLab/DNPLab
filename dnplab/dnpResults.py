import matplotlib.pyplot as plt
import numpy as np

from .dnpData import dnpdata

figure = plt.figure
legend = plt.legend
xlim = plt.xlim
ylim = plt.ylim
gca = plt.gca

dark_green = "#46812B"
light_green = "#67AE3E"
dark_grey = "#4D4D4F"
light_grey = "#A7A9AC"
orange = "#F37021"

plt.rcParams["lines.linewidth"] = 1.5
plt.rcParams["axes.prop_cycle"] = plt.cycler(
    color=[orange, dark_green, light_green, dark_grey, light_grey]
)


def imshow(data, *args, **kwargs):  # TODO: drop unused args and kwargs
    """Image Plot for dnpdata object

    Args:
        data (dnpdata): dnpdata object for image plot
        args: args for matplotlib imshow function
        kwargs: kwargs for matplotlib imshow function

    Example::

       # Plotting a dnpdata object
       dnp.dnpResults.plt.figure()
       dnp.dnpResults.imshow(data)
       dnp.dnpResults.plt.show()

       # Plotting a workspace (dnpdata_collection)
       dnp.dnpResults.plt.figure()
       dnp.dnpResults.imshow(ws['proc'])
       dnp.dnpResults.plt.show()
    """

    dims = data.dims

    x_coord = data.coords[dims[1]]
    y_coord = data.coords[dims[0]]

    x_min = np.min(x_coord)
    x_max = np.max(x_coord)
    y_min = np.min(y_coord)
    y_max = np.max(y_coord)

    plt.imshow(data.values, aspect="auto", extent=[x_min, x_max, y_max, y_min])
    plt.xlabel(dims[1])
    plt.ylabel(dims[0])


def stack(data, *args, offset=None, **kwargs):
    """Stack Plot for 2D data

    Args:
        data (dnpdata): dnpdata object for matplotlib plot function
        args: args for matplotlib plot function
        offset: Value to offset each spectra, by default maximum of absolute value
        kwargs: kwargs for matplotlib plot function

    Example::

       dnp.dnpResults.plt.figure()
       dnp.dnpResults.stack(data)
       dnp.dnpResults.plt.show()

    """

    coord = data.coords[0]
    dim = data.dims[0]

    if offset == None:
        offset = np.max(data.abs)

    offset_matrix = (
        offset * np.ones(coord.size).reshape(-1, 1) * np.r_[0 : data.coords[1].size]
    )

    plt.plot(coord, data.values + offset_matrix, *args, **kwargs)


def waterfall(data, dx, dy, *args, **kwargs):
    """Waterfall plot for 2d data

    Args:
        data (dnpData): 2d Data object for waterfall plot
        dx (float, int): x-increment for each line
        dy (float, int): y-increment for each line

    Example::

       dnp.dnpResults.plt.figure()
       dnp.dnpResults.waterfall(data)
       dnp.dnpResults.plt.show()

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


def plot(data, dim=None, *args, **kwargs):
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
    if not dim:
        dim = data.dims[0]

    coord = data.coords[dim]

    plt.plot(coord, data.values.real, *args, **kwargs)


show = plt.show
