import numpy as np

from dnplab import dnpData


def save_csv(filename, odnpData: dnpData):
    data = odnpData.values
    np.savetxt(filename, data, delimiter=",")
