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


show = plt.show