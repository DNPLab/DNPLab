import numpy as np
import os


def load_DAT(filename):
    '''Import Bruker Binary .DAT ELEXSYS data

    Args:
        filename (str): path to data file in string format

    Returns:
        numpy.ndarray: Experimental data as numpy array. No reshaping or formatting applied. Complex or multi-dimensional data must be subsequently processed.
    '''
    
    data = np.fromfile(filename, dtype = '>f8')

    return data

def load_DSC(filename):
    '''Import Bruker ASCII Parameters .DSC File

    Args:
        filename (str): Path to parameters file in string format

    Returns:
        dict: Dictionary of Descriptor Information in DSC file
    '''

    DSC_dict = {}
    with open(filename,'r') as f:
        line = f.readline()

        if line[0:5] == '#DESC':
            desc_line = '*'
            while desc_line[0] != '#':
                desc_line = f.readline()
                if (desc_line[0] != '#') and (desc_line[0] != '*'):
                    key, value = desc_line.rstrip().rsplit('\t')

                    # Remove apostrophes surrounding certain variables
                    if (value[0] == "'") and (value[-1] == "'"):
                        value = value[1:-1]

                    # Attempt to convert to float or int
                    try:
                        if '.' in value:
                            value = float(value)
                        else:
                            value = int(value)
                    except:
                        pass

                    DSC_dict[key] = value

    return DSC_dict

def load_elexsys(filename):
    '''Import Elexsys Data

    Args:
        filename (str): Filename of data

    Returns: 
        tuple: tuple containing:

            dimsList (*numpy.ndarray*, *list*): Numpy array of data if data is 1d, otherwise returns a list of the dimension arrays

            data (*numpy.ndarray*): Numpy array of data
    '''
    # Handle extension
    filename,extension = os.path.splitext(filename)

    raw_data = load_DAT(filename + '.DTA')

    DSC_dict = load_DSC(filename + '.DSC')

    if DSC_dict['IKKF'] == 'CPLX':
        data = raw_data[0::2] + 1j * raw_data[1::2]


    # Set points equal to one for reshaping
    xpts = 1
    ypts = 1
    zpts = 1

    dimsList = []
    dims = 0
    if DSC_dict['XTYP'] != 'NODATA':
        xpts = DSC_dict['XPTS']
        xmin = DSC_dict['XMIN']
        xwid = DSC_dict['XWID']

        x = xmin + np.r_[0:xwid:1j*xpts]

        dims += 1
        dimsList.append(x)

    if DSC_dict['YTYP'] != 'NODATA':
        ypts = DSC_dict['YPTS']
        ymin = DSC_dict['YMIN']
        ywid = DSC_dict['YWID']

        y = ymin + np.r_[0:ywid:1j*ypts]
        dims += 1
        dimsList.append(y)

    if DSC_dict['ZTYP'] != 'NODATA':
        zpts = DSC_dict['ZPTS']
        zmin = DSC_dict['ZMIN']
        zwid = DSC_dict['ZWID']

        z = zmin + np.r_[0:zwid:1j*zpts]
        dims += 1
        dimsList.append(z)

    # Reshape data
    data = data.reshape(xpts,ypts,zpts)

    # Remove length 1 dimensions
    data = data.squeeze() 
    if len(dimsList) == 1:
        dimsList = dimsList[0]

    return dimsList, data

if __name__ == '__main__':
    pass
