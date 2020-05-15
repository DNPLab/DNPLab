from .. import odnpData
import numpy as np

def save_csv(odnpData,filename):
    '''
    '''

    data = odnpData.data
    np.savetxt(data,delimiter = ',')


