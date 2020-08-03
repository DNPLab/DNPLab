import matplotlib.pyplot as plt
import numpy as np

#import dnpData
from . import dnpData
#import dnpData

figure = plt.figure
legend = plt.legend
xlim = plt.xlim
ylim = plt.ylim
gca = plt.gca

dark_green = '#46812B'
light_green = '#67AE3E'
dark_grey = '#4D4D4F'
light_grey = '#A7A9AC'
orange = '#F37021'

plt.rcParams['lines.linewidth'] = 1.5
plt.rcParams['axes.prop_cycle'] = plt.cycler(color = [orange, dark_green, light_green, dark_grey, light_grey])

def plot(data, *args, **kwargs):
    '''
    '''
#    coord = data.coords[dim]
    coord = data.coords[0]
    dim = data.dims[0]

#    original_order = data.dims
#    data.reorder([dim])
    plt.plot(coord, data.values, *args, **kwargs)
    plt.xlabel(dim)
#    data.reorder(original_order)


def imshow(data, *args, **kwargs):
    '''
    '''

    dims = data.dims

    x_coord = data.coords[dims[1]]
    y_coord = data.coords[dims[0]]

    x_min = np.min(x_coord)
    x_max = np.max(x_coord)
    y_min = np.min(y_coord)
    y_max = np.max(y_coord)

    plt.imshow(data.values, aspect = 'auto', extent = [x_min, x_max, y_max, y_min])
    plt.xlabel(dims[1])
    plt.ylabel(dims[0])


#    return NotImplemented

show = plt.show

if __name__ == '__main__':

    x = np.r_[-10:10:100j].reshape(10,10)
    y = x**2.

    data = dnpData.dnpData(y, [np.r_[0:10], np.r_[0:10]], ['x'])
    plot(data)
    show()

#    fig = plt.figure()
#    plt.plot(x, y)
#
