import numpy as np


def import_winepr(filename):
    '''
    '''
    return

def load_par(filename):
    '''Import WinEPR ".par" parameters file

    Args:
        filename (str): Path to .par file

    Returns:
        dict: Dictionary of parameters
    '''

    with open(filename, 'r') as f:
        raw = f.read()

    lines = raw.rstrip().rsplit('\n')

    params_dict = {}

    for line in lines:
        split_line = line.split(' ', 1)
        key = split_line[0]

        value = split_line[1]
        value = value.replace(' ', '')

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

        params_dict[key] = value

    return params_dict

def load_spc(filename):
    '''Import WinEPR Binary .spc file

    Args:
        filename (str): Path to file

    Returns:
        numpy.ndarray
    '''

    data = np.fromfile(filename, '<f')

    return data

def load_winepr(filename):
    '''Import WinEPR data

    Args:
        filename (str): Filename of data

    Returns:
        tuple: tuple containing:

            dimsList (*numpy.ndarray*, *list*): Numpy array of data if data is 1d, otherwise returns a list of the dimension arrays

            data (*numpy.ndarray*): Numy array of data

    '''

    # Handle extension

    filename, extension = os.path.splitext(filename)

    data = load_spc(filename + '.spc')

    params_dict = load_par(filename + '.par')

    center_field = params_dict['HCF']
    field_width = params_dict['HSW']
    pts = params_dict['ANZ']

    x = np.r[-1*field_width/2:field_width/2:1j*pts] + center_field

    return x, data

