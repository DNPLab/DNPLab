from .. import odnpData
import numpy as np
from struct import unpack
import os
import glob


def importKea(path, parameters_filename = None, verbose = False):
    '''Import Kea data

    Args:
        path (str): Path to data
        num (int): Experiment number
        verbose (bool): If true, prints additional information for troubleshooting
    
    Returns:
        odnpData object with Kea data
    '''

    if parameters_filename == None:
        parameters_filename = 'acqu.par'
    extension = ''

    if os.path.isfile(path):
        path, filename = os.path.split(path)
        filename, extension = os.path.splitext(filename)
    elif os.path.isdir(path):
        filesList = glob.glob(path + '*.[1-4]d')
        if len(filesList) == 0:
            raise ValueError('No binary data file in directory:')
        elif len(filesList) > 1:
            raise ValueError('More than one binary data file in directory:',fileList)
        else:
            data_filename = filesList[0]

    attrs = import_par(os.path.join(path, parameters_filename))

    if extension == '.csv':
        # Import csv data
        x, data = import_csv(os.path.join(path, filename + extension))
    else:
        # Import Binary data
        x, data = import_nd(os.path.join(path, filename + extension))

    if 'b1Freq' in attrs:
        nmrFrequency = attrs['b1Freq']
    else:
        nmrFrequency = 14.5
    try:
        nmrFrequency = float(nmrFrequency)
    except:
        nmrFrequency = float(nmrFrequency.replace('d',''))

    attrs['nmrFrequency'] = nmrFrequency

    # Assume direct dimension is 1st dimension
    data_shape = np.shape(np.squeeze(data))

    dims = []
    coords = []
    for ix in range(len(data_shape)):
        dims.append(str(ix))
        coords.append(np.arange(data_shape[ix]))

    # If axes information is give, assume it is the first dimension
    if x is not None:
        dims[0] = 't'
        coords[0] = x

    kea_data = odnpData(data, coords, dims, attrs)

    return kea_data

def import_nd(path):
    '''Import Kea 1d, 2d, 3d, 4d files

    Args:
        path (str): Path to file

    Returns:
        tuple:
            x (None, numpy.array): Axes if included in binary file, None otherwise
            data (numpy.array): Numpy array of data
    '''

    ascii_header_size = 12
    ascii_header_format = '<4s4s4s'

    with open(path, 'rb') as f:
        ascii_header_bytes = f.read(ascii_header_size)

        ascii_header = unpack(ascii_header_format, ascii_header_bytes)

        owner = ascii_header[0][::-1].decode('utf-8')
        format_ = ascii_header[1][::-1].decode('utf-8')
        version = ascii_header[2][::-1].decode('utf-8')

        # header size depends on version
        if version == 'V1.0':
            header_format = '4i'
            header_size = 16
        else:
            header_format = '5i'
            header_size = 20

        header_bytes = f.read(header_size)

        header = unpack(header_format, header_bytes)

        dataType = header[0]
        xDim = header[1]
        yDim = header[2]
        zDim = header[3]
        if version != 'V1.0':
            qDim = header[4]
        else:
            qDim = 1

        raw = f.read()

        x = None
        if dataType == 500: # float
            raw_data = unpack('<%if'%(xDim*yDim*zDim*qDim), raw)
            data = np.array(raw_data)
        elif dataType == 501: # complex
            raw_data = unpack('<%if'%(xDim*yDim*zDim*qDim*2), raw)
            data = np.array(raw_data)
            data = data[0::2] + 1j * data[1::2]
        elif dataType == 502: # double
            raw_data = unpack('<%id'%(xDim*yDim*zDim*qDim), raw)
            data = np.array(raw_data)
        elif dataType == 503:
            raw_data = unpack('<%if'%(xDim*yDim*zDim*qDim*2), raw)
            raw_data = np.array(raw_data)
            x = raw_data[0:xDim]
            data = raw_data[xDim:]
        elif dataType == 504:
            raw_data = unpack('<%if'%(xDim*yDim*zDim*qDim*3), raw) #504
            raw_data = np.array(raw_data)
            x = raw_data[0:xDim]
            data = raw_data[xDim:]
            data = data[0::2] + 1j * data[1::2]
        else:
            raise ValueError('Data %i type not recognized'%dataType)

        data = data.reshape(xDim, yDim, zDim, qDim) # reshape data
        data = data.squeeze() # remove length 1 dimensions

    return x, data

def import_par(path):
    ''' Import Kea parameters .par file

    Args:
        path (str): Path to parameters file

    Returns:
        dict: Dictionary of Kea Parameters
    '''

    attrs = {}

    with open(path, 'r') as f:

        raw = f.read()

        lines = raw.rstrip().rsplit('\n')

        for line in lines:
            key, value = line.rstrip().rsplit(' = ')

            if value[0] == '"' and value[-1] == '"':
                value = value[1:-1]

            if '.' in value:
                try:
                    value = float(value)
                except:
                    pass
            else:
                try:
                    value = int(value)
                except:
                    pass
            attrs[key] = value

    return attrs 

def import_csv(path, return_raw = False, is_complex = True):
    '''Import Kea csv file

    Args:
        path (str): Path to csv file

    Returns:
        tuple:
            x(numpy.array): axes if return_raw = False
            data(numpy.array): Data in csv file
    '''

    raw = np.loadtxt(path, delimiter = ',')

    if not return_raw:
        x = raw[:,0]
        if is_complex:
            data = raw[:,1::2] + 1j * raw[:,2::2]
        else:
            data = raw[:,1:]
        data = np.squeeze(data)
        return x, data
    else:
        return raw
