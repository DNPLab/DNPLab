from .. import dnpData
import numpy as np
import h5py


def saveh5(dataDict, path, overwrite = False):
    '''
    Save All Data in .h5 format
    Requires a dictionary of dnpData objects

    '''

    if overwrite:
        mode = 'w'
    else:
        mode = 'w-'

    keysList = dataDict.keys()

    f = h5py.File(path, mode)

#    f.attrs['dnpLab_version'] = version

    for key in keysList:
        dnpDataObject = dataDict[key]
        
        dnpDataGroup = f.create_group(key,track_order = True)
        dnpDataGroup.attrs['dnpLab_version'] = dataDict[key].version
        dims_group = dnpDataGroup.create_group('dims') # dimension names e.g. x,y,z
        attrs_group = dnpDataGroup.create_group('attrs') # dictionary information
        dnp_dataset = dnpDataGroup.create_dataset('values',data = dnpDataObject.data)

        # Save axes information
        for ix in range(len(dnpDataObject.axes)):
            label = dnpDataObject.axesLabels[ix]
            this_axes = dnpDataObject.axes[ix]
            dims_group.create_dataset(label,data = this_axes)
            dims_group[label].make_scale(label)

            dnp_dataset.dims[ix].attach_scale(dims_group[label])

        # Save Parameters
        for key in dnpDataObject.params:
            attrs_group.attrs[key] = dnpDataObject.params[key]
    f.close()

def loadh5(path):
    '''
    Returns Dictionary of dnpDataObjects
    '''

    dnpDict = {}

    f = h5py.File(path,'r')
    keysList = f.keys()
    print('keys:')
    print(keysList)

    for key in keysList:
        axes = []
        axesLabels = []
        params = {}
        data = f[key]['values'][:]
        version = f[key].attrs['dnpLab_version']

        for index in range(len(np.shape(data))):
            dimKey = f[key]['values'].dims[index].keys()[0] # assumes 1 key only
            axes.append(f[key]['values'].dims[index][dimKey][:])
            axesLabels.append(dimKey)

        for k in f[key]['attrs'].attrs.keys():
            print(k)
            print(f[key]['attrs'].attrs[k])
            params[k] = f[key]['attrs'].attrs[k]
        dnpDict[key] = dnpData(data,axes,axesLabels,params)
        dnpDict[key].version = version

    return dnpDict

