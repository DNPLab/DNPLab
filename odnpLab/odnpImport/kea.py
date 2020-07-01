from .. import odnpData
import numpy as np
from struct import unpack


def importKea(path,filename = '',num = 1 ,verbose = False):
    '''Import Kea data

    Args:
        path (str): Directory of data
        filename (str): Filename of data
        num (int): Experiment number
        verbose (bool): If true, prints additional information for troubleshooting
    
    Returns:
        odnpData object with Kea data
    '''
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
    '''Load EMX data
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

def importd01(filename):
    '''Import Prospa Binary Format

    Args:
        filename: Path to data

    Returns:
        numpy.ndarray: binary data
    '''

    headerSize = 32

    header_fmt = '<8i'
    with open(filename,'rb') as f:
        headerString = f.read(headerSize)
        header = unpack(header_fmt,headerString)

        data_amount = header[0]
        data_format = header[1]
        data_dimension = header[2]
        array_of_dimensions = header[3:7]
        overall_data_size = header[7]

        if data_format == 0:
            unpack_type = 'd'
            bytes_per_point = 8
        else:
            unpack_type = 'f'
            bytes_per_point = 4



        total_bytes = bytes_per_point*overall_data_size
        dataString = f.read(4)

        data = unpack('<f',dataString)

    return data


