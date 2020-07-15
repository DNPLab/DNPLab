from .. import dnpData
import numpy as np

def save_csv(filename,odnpData):
    '''
    '''

    data = odnpData.data

    np.savetxt(filename,data,delimiter = ',')


