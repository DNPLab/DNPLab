import matplotlib.pyplot as plt
import numpy as np

#import dnpData
from . import dnpData
#import dnpData

figure = plt.figure

def plot(data, dim, *args, **kwargs):
    '''
    '''
    coord = data.coords[dim]
    original_order = data.dims
    data.reorder([dim])
    plt.plot(coord, data.values, *args, **kwargs)
    plt.xlabel(dim)
    data.reorder(original_order)


def imshow(data, *args, **kwargs):
    '''
    '''

    dims = data.dims

    x_coord = data.coords[dims[0]]
    y_coord = data.coords[dims[1]]

    x_min = np.min(x_coord)
    x_max = np.max(x_coord)
    y_min = np.min(y_coord)
    y_max = np.max(y_coord)

    plt.imshow(data.values, aspect = 'auto', extent = [x_min, x_max, y_min, y_max])
    plt.xlabel(dims[0])
    plt.ylabel(dims[1])


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
