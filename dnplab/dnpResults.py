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
    coord = data.coords[0]
    dim = data.dims[0]

    plt.plot(coord, data.values, *args, **kwargs)


show = plt.show
