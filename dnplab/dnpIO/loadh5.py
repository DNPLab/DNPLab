from .. import dnpdata, create_workspace
import warnings
import numpy as np
import h5py


def load_h5(path):
    """Returns Dictionary of dnpDataObjects

    Args:
        path (str): Path to h5 file

    Returns:
        dnpdata_collection: workspace object with data
    """

    dnp_dict = {}

    f = h5py.File(path, "r")
    keys_list = f.keys()

    for key in keys_list:

        if f[key].attrs["dnplab_data_type"] == "dnpdata":
            data = read_dnpdata(f[key])
        elif f[key].attrs["dnplab_data_type"] == "dict":
            data = read_dict(f[key])
        else:
            warnings.warn("could not import key: %s" % str(key))

        dnp_dict[key] = data
    return create_workspace(dnp_dict)


def read_dnpdata(dnpdata_group):
    axes = []
    dims = []
    attrs = {}
    values = dnpdata_group["values"][:]
    version = dnpdata_group.attrs["dnplab_version"]

    for index in range(len(np.shape(values))):
        dim_key = dnpdata_group["values"].dims[index].keys()[0]  # assumes 1 key only
        axes.append(dnpdata_group["values"].dims[index][dim_key][:])
        dims.append(dim_key)

    for k in dnpdata_group["attrs"].attrs.keys():
        attrs[k] = dnpdata_group["attrs"].attrs[k]

    data = dnpdata(values, axes, dims, attrs)

    if "proc_attrs" in dnpdata_group.keys():
        proc_attrs = []
        for k in dnpdata_group["proc_attrs"].keys():
            proc_attrs_name = k.split(":", 1)[1]
            proc_attrs_dict = dict(dnpdata_group["proc_attrs"][k].attrs)
            data.add_proc_attrs(proc_attrs_name, proc_attrs_dict)
    return data


def read_dict(dnpdata_group):
    data = dict(dnpdata_group["attrs"].attrs)
    return data
