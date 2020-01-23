from .. import odnpData
import numpy as np
import h5py

def saveh5(dataDict,path):
    '''
    Save All Data in .h5 format
    Requires a dictionary of odnpData objects

    '''

    if not isinstance(dataDict,dict):
        print('data input must be dictionary of odnpData objects')
        return

    keysList = dataDict.keys()

    f = h5py.File(path,'w')
    for key in keysList:
        odnpDataObject = dataDict[key]

        h5_data = f.create_dataset(key,data = odnpDataObject.data)
        for index in range(len(odnpDataObject.axesLabels)):
            label = odnpDataObject.axesLabels[index]
            detailedLabel = '*' + key + ' ' + label + ' axes*'
            axes = odnpDataObject.axes[index].copy()

            f[detailedLabel] = axes
            f[detailedLabel].make_scale(detailedLabel)
            f[key].dims[index].attach_scale(f[detailedLabel])

            for k in odnpDataObject.params.keys():
                f[key].attrs[k] = odnpDataObject.params[k]
    f.close()

def loadh5(path):
    '''
    Returns Dictionary of odnpDataObjects
    '''

    odnpDict = {}

    f = h5py.File(path,'r')

    keysList = f.keys()
    dataKeysList = []
    for key in keysList:
        if (key[-1] == '*') and (key[0] == '*'):
            pass
            # disregard as data, assume axes
        else:
            dataKeysList.append(key)
    
    for key in dataKeysList:
        data = f[key][:]
        params =  {} 
        for k in f[key].attrs.keys():
            if k != 'DIMENSION_LIST':
                params[k] = f[key].attrs[k]

        axesLabels = []
        axes = []

        for index in range(len(np.shape(data))):
            dimKeyRaw = f[key].dims[index].keys()[0] # assumes 1 key only
            label = dimKeyRaw[len('*' + key + ' '):-1*len(' axes*')]

            axesLabels.append(label)
            axes.append(f[key].dims[index][dimKeyRaw][:])
        odnpDict[key] = odnpData(data,axes,axesLabels,params)
    

    f.close()
    return odnpDict



