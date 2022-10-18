import matplotlib.pyplot as plt
import numpy as np


def imshow(data, *args, **kwargs):  # TODO: drop unused args and kwargs
    """Image Plot for dnpdata object

    Args:
        data (DNPData): DNPData object for image plot
        args: args for matplotlib imshow function
        kwargs: kwargs for matplotlib imshow function

    Returns:
        Returns formated matplotlib plot.

    Example:

       Plotting a dnpdata object

       >>> dnp.plt.figure()
       >>> dnp.imshow(data)
       >>> dnp.plt.show()

       Plotting a workspace (dnpdata_collection)

       >>> dnp.plt.figure()
       >>> dnp.imshow(data)
       >>> dnp.plt.show()
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
