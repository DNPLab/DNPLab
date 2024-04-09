import os
import re
import numpy as _np
from matplotlib.pylab import *
from ..core.util import concat
from .. import DNPData

try:
    import skrf as _rf

except Exception as e:
    raise ImportError("Please install scikit-rf first.")


def import_vna(path, *args, **kwargs):
    """Import VNA data and return dnpdata object"""
    x, values, attrs, dim = import_snp(path, *args, **kwargs)
    dnpDataObject = get_dnpdata(values, x, attrs, dim)
    return dnpDataObject


def get_dnpdata(values, coords, attrs, concat_dim=None):
    if len(_np.shape(values)) == 1:
        return DNPData(
            values, coords=[coords], dims=["f"], attrs=attrs, dnplab_attrs=attrs
        )

    else:
        new_coords = list(range(_np.shape(values)[0]))
        data_list = []
        for value in values:
            data = get_dnpdata(value, coords=coords, attrs=attrs)
            data_list.append(data)

        return concat(data_list, concat_dim, new_coords)


def import_snp(path, *args, **kwargs):
    """Import sNp file and return numpy array"""
    _, extension = os.path.splitext(path)

    extension_reg_ex = "[.]s[0-9]{1,}p"
    if re.fullmatch(extension_reg_ex, extension) == None:
        raise ValueError("File Extension Not Given, Unspecified sNp file")

    num_reg_ex = "[0-9]{1,}"
    num = int(re.search(num_reg_ex, extension)[0])

    if num > 2:
        raise ValueError("Currently on s1p and s2p are supported")

    data = _rf.Network(path, *args, **kwargs)

    attrs = {
        "data_format": "VNA",
        "center_field": data.frequency.center,
        "sweep_field": data.frequency.span,
    }

    if extension == ".s1p":
        value = data.s[:, 0, 0]
        x = data.f
        dim = None
        attrs["data_order"] = ["s11"]

    elif extension == ".s2p":
        s11 = data.s[:, 0, 0]
        s12 = data.s[:, 0, 1]
        s21 = data.s[:, 1, 0]
        s22 = data.s[:, 1, 1]
        value = _np.array([s11, s12, s21, s22])
        x = data.f
        dim = "s"
        attrs["data_order"] = ["s11", "s12", "s21", "s22"]

    return x, value, attrs, dim
