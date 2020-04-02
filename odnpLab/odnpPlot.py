from . import odnpData as _odnpData

import numpy as _np

import matplotlib.pylab as plt

def plot(odnpData,axesLabel,*args,**kwargs):
    '''Plot data with given dimension as X-axes
    '''

#    if 'axesLabel' in *args:
#        axesLabel == *args['axesLabel']
#        *args.pop('axesLabel')
#
    index = odnpData.index(axesLabel)

    data = odnpData.data

    plt.plot(odnpData.axes[index],_np.swapaxes(odnpData.data,0,index),*args,**kwargs)
    plt.xlabel(odnpData.axesLabels[index])
