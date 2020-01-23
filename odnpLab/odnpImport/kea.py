#import odnpData
#from odnpData import odnpData
from .. import odnpData
import numpy as np


def importKea(path,filename = '',num = 1 ,verbose = False):
    params_dict = {}
    try:
        with open(path + filename + '/%i/'%num + 'acqu.par','r') as f:
            raw_params = f.read()


        raw_params = raw_params.strip().split('\n')
        for line in raw_params:
            if verbose:
                print(line)
            key_value = line.split(' = ')
            try:
                params_dict[key_value[0]] = float(key_value[1])
            except:
                params_dict[key_value[0]] = key_value[1]
        nmr_params = {}
        nmrFreq = float(params_dict['b1Freq'].strip('d'))*1e6 # convert to Hz

        params_dict['nmrFreq'] = nmrFreq
    except:
        pass

    raw_data = np.loadtxt(path + filename + '/%i/'%num + 'data.csv', delimiter = ',')

    t = raw_data[:,0]
    t = t / 1.e6 # convert from us to s

    temp_data = raw_data[:,1] - 1j*raw_data[:,2]
#    data = odnpData.odnpData(temp_data,[t],['t'],nmr_params)
#    data = odnpData(temp_data,[t],['t'],nmr_params)
    data = odnpData(temp_data,[t],['t'],params_dict)
#    data.kea_params = params_dict
    return data

def importEmx(path,filename = ''):
    '''
    Load EMX data
    '''

    if path[-1] != '\\' and path[-1] != '/':
        filenamePath = path + filename
    else:
        filenamePath = path + '/' + filename

    # Strip Extension from filename path string
    if (filenamePath[-4:] == '.spc') or (filenamePath[-4:] == '.par'):
        ext = filenamePath[-4:]
        filenamePath = filenamePath[:-4]

    data = np.fromfile(filenamePath + '.spc',dtype = np.float32)

    # open with 'U' mode to enable univeral new-line for Python2 only
    with open(filenamePath + '.par',mode = 'U') as f:
        params_raw = f.read()
    params_lines = params_raw.strip().split('\n')

    paramsDict = {}

    for line in params_lines:
        key = line[0:3]
        value = line[4:]
#        print(key, value)

        try:
            value = float(value)
        except:
            pass

        paramsDict[key] = value
    print(paramsDict)

    center_field = paramsDict['HCF']
    sweep_width = paramsDict['HSW']
    xpts = paramsDict['ANZ']


    field = center_field + np.r_[-0.5*sweep_width:0.5*sweep_width:1j*xpts]

    output = odnpData(data,[field],['field'],paramsDict)
    return output

#def saveh5(dataDict,path):
#    '''
#    Save All Data in .h5 format
#    Requires a dictionary of odnpData objects
#
#    '''
#
#    if not isinstance(dataDict,dict):
#        print('data input must be dictionary of odnpData objects')
#        return
#
#    keysList = dataDict.keys()
#
#    f = h5py.File(path,'w')
#    for key in keysList:
#        odnpDataObject = dataDict[key]
#
#        h5_data = f.create_dataset(key,data = odnpDataObject.data)
#        for index in range(len(odnpDataObject.axesLabels)):
#            label = odnpDataObject.axesLabels[index]
#            detailedLabel = '*' + key + ' ' + label + ' axes*'
#            axes = odnpDataObject.axes[index].copy()
#
#            f[detailedLabel] = axes
#            f[detailedLabel].make_scale(detailedLabel)
#            f[key].dims[index].attach_scale(f[detailedLabel])
#
#            for k in odnpDataObject.params.keys():
#                f[key].attrs[k] = odnpDataObject.params[k]
#    f.close()
#
#def loadh5(path):
#    '''
#    Returns Dictionary of odnpDataObjects
#    '''
#
#    odnpDict = {}
#
#    f = h5py.File(path,'r')
#
#    keysList = f.keys()
#    dataKeysList = []
#    for key in keysList:
#        if (key[-1] == '*') and (key[0] == '*'):
#            pass
#            # disregard as data, assume axes
#        else:
#            dataKeysList.append(key)
#    
#    for key in dataKeysList:
#        data = f[key][:]
#        params =  {} 
#        for k in f[key].attrs.keys():
#            if k != 'DIMENSION_LIST':
#                params[k] = f[key].attrs[k]
#
#        axesLabels = []
#        axes = []
#
#        for index in range(len(np.shape(data))):
#            dimKeyRaw = f[key].dims[index].keys()[0] # assumes 1 key only
#            label = dimKeyRaw[len('*' + key + ' '):-1*len(' axes*')]
#
#            axesLabels.append(label)
#            axes.append(f[key].dims[index][dimKeyRaw][:])
#        odnpDict[key] = odnpData(data,axes,axesLabels,params)
#    
#
#    f.close()
#    return odnpDict
#
#
#
