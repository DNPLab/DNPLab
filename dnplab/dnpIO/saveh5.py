from .. import dnpdata
import warnings
import numpy as np
import h5py


def save_h5(dataDict, path, overwrite=False):
    """Save workspace in .h5 format

    Args:
        dataDict (dnpdata_collection): dnpdata_collection object to save.
        path (str): Path to save data
        overwrite (bool): If True, h5 file can be overwritten. Otherwise, h5 file cannot be overwritten
    """

    if overwrite:
        mode = "w"
    else:
        mode = "w-"

    keysList = dataDict.keys()

    f = h5py.File(path, mode)

    for key in keysList:
        dnpDataObject = dataDict[key]

        dnpDataGroup = f.create_group(key, track_order=True)
        if isinstance(dnpDataObject, dnpdata):
            write_dnpdata(dnpDataGroup, dnpDataObject)
        elif isinstance(dnpDataObject, dict):
            write_dict(dnpDataGroup, dnpDataObject)
        else:
            warnings.warn("Could not write key: %s" % str(key))

    f.close()


def write_dnpdata(dnpDataGroup, dnpDataObject):
    """Takes file/group and writes dnpData object to it

    Args:
        dnpDataGroup: h5 group to save data to
        dnpDataObject: dnpdata object to save in h5 format
    """
    dnpDataGroup.attrs["dnplab_version"] = dnpDataObject.version
    dnpDataGroup.attrs["dnplab_data_type"] = "dnpdata"
    dims_group = dnpDataGroup.create_group("dims")  # dimension names e.g. x,y,z
    attrs_group = dnpDataGroup.create_group("attrs")  # dictionary information
    dnp_dataset = dnpDataGroup.create_dataset("values", data=dnpDataObject.values)

    # Save axes information
    for ix in range(len(dnpDataObject.coords)):
        label = dnpDataObject.dims[ix]
        this_axes = dnpDataObject.coords[ix]
        dims_group.create_dataset(label, data=this_axes)
        dims_group[label].make_scale(label)

        dnp_dataset.dims[ix].attach_scale(dims_group[label])

    # Save Parameters
    for key in dnpDataObject.attrs:
        attrs_group.attrs[key] = dnpDataObject.attrs[key]

    # Save proc_steps
    if hasattr(dnpDataObject, "proc_attrs"):
        proc_attrs = dnpDataObject.proc_attrs
        proc_attrs_group = dnpDataGroup.create_group("proc_attrs", track_order=True)
        for ix in range(len(proc_attrs)):
            proc_step_name = proc_attrs[ix][0]
            proc_dict = proc_attrs[ix][1]
            proc_attrs_group_subgroup = proc_attrs_group.create_group(
                "%i:%s" % (ix, proc_step_name)
            )
            for key in proc_dict:
                value = proc_dict[key]
                proc_attrs_group_subgroup.attrs[key] = value


def write_dict(dnpDataGroup, dnpDataObject):
    """Writes dictionary to h5 file"""
    #    dnpDataGroup.attrs['dnplab_version'] = dnpDataObject.version
    dnpDataGroup.attrs["dnplab_data_type"] = "dict"
    #    dnpDataGroup.attrs['dnplab_version'] = dnpDataObject.version
    attrs_group = dnpDataGroup.create_group("attrs")

    for key in dnpDataObject.keys():
        attrs_group.attrs[key] = dnpDataObject[key]
