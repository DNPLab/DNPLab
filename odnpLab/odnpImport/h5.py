from .. import odnpData
import numpy as np
import h5py

odnp_data_h5version = '1.0'

def saveh5(dataDict, path, overwrite = False):
    '''
    Save All Data in .h5 format
    Requires a dictionary of odnpData objects

    '''

    if overwrite:
        mode = 'w'
    else:
        mode = 'w-'

    keysList = dataDict.keys()

    f = h5py.File(path, mode)

    for key in keysList:
        odnpDataObject = dataDict[key]
        odnpDataGroup = f.create_group(key,track_order = True)
        dims_group = odnpDataGroup.create_group('dims') # dimension names e.g. x,y,z
        attrs_group = odnpDataGroup.create_group('attrs') # dictionary information
        odnp_dataset = odnpDataGroup.create_dataset('values',data = odnpDataObject.data)

        # Save axes information
        for ix in range(len(odnpDataObject.axes)):
            label = odnpDataObject.axesLabels[ix]
            this_axes = odnpDataObject.axes[ix]
            dims_group.create_dataset(label,data = this_axes)
            dims_group[label].make_scale(label)

            odnp_dataset.dims[ix].attach_scale(dims_group[label])

        # Save Parameters
        for key in odnpDataObject.params:
            attrs_group.attrs[key] = odnpDataObject.params[key]
    f.close()

def loadh5(path):
    '''
    Returns Dictionary of odnpDataObjects
    '''

    odnpDict = {}

    f = h5py.File(path,'r')
    keysList = f.keys()
    print('keys:')
    print(keysList)

    for key in keysList:
        axes = []
        axesLabels = []
        params = {}
        data = f[key]['values'][:]

        for index in range(len(np.shape(data))):
            dimKey = f[key]['values'].dims[index].keys()[0] # assumes 1 key only
            axes.append(f[key]['values'].dims[index][dimKey][:])
            axesLabels.append(dimKey)

        for k in f[key]['attrs'].attrs.keys():
            print(k)
            print(f[key]['attrs'].attrs[k])
            params[k] = f[key]['attrs'].attrs[k]
        odnpDict[key] = odnpData(data,axes,axesLabels,params)

    return odnpDict

