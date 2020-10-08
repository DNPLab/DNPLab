# TODO: remove unused imports
import numpy as np
import os
import re
from matplotlib.pylab import *

from .. import dnpdata as _dnpdata


def import_vna(path):
    """Import VNA data and return dnpdata object"""
    x, data = import_snp(path)
    # Not General
    dnpDataObject = _dnpdata(data, [x], ["f"], {})

    return dnpDataObject


# TODO: remove prints or make them optional
def import_snp(path):
    """Import sNp file and return numpy array"""
    path_filename, extension = os.path.splitext(path)

    extension_reg_ex = "[.]s[0-9]{1,}p"
    print(re.fullmatch(extension_reg_ex, extension))
    print(extension)
    if re.fullmatch(extension_reg_ex, extension) == None:
        raise ValueError("File Extension Not Given, Unspecified sNp file")

    num_reg_ex = "[0-9]{1,}"
    num = int(re.search(num_reg_ex, extension)[0])

    print(num)
    if num > 2:
        raise ValueError("Currently on s1p and s2p files are supported")

    f = open(path)
    read_string = " "
    while read_string[0] != "#":
        read_string = f.readline()
    raw = np.genfromtxt(f, skip_header=2, defaultfmt="11f")
    f.close()

    if num == 1:
        x = raw[:, 0]
        data = raw[:, 1] + 1j * raw[:, 2]

    if num == 2:
        x = raw[:, 1]

        data = np.zeros((len(x), 2, 2))

        data[:, 0, 0] = raw[:, 1] + 1j * raw[:, 2]  # S11
        data[:, 1, 0] = raw[:, 3] + 1j * raw[:, 4]  # S21
        data[:, 0, 1] = raw[:, 5] + 1j * raw[:, 6]  # S12
        data[:, 1, 1] = raw[:, 7] + 1j * raw[:, 8]  # S22

    if num > 2:
        x = raw[0::num]
        data = np.zeros((len(x), num, num))

        # TODO: Use list comprehension instead of two for loops
        for n in range(num):

            for m in range(num):

                data[:, n, m] = raw[n::num, 1 + 2 * m] + 1j * raw[n::num, 2 * (1 + m)]

    return x, data
