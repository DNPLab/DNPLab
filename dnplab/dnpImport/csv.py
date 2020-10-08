from .. import dnpData
import numpy as np


def save_csv(filename, odnpData: dnpData):
    data = odnpData.values
    np.savetxt(filename, data, delimiter=",")
