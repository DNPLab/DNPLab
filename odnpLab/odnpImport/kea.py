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


def import_nd(path):
    '''Import Kea 1d, 2d, 3d files

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
        return x, data
    else:
        return raw
