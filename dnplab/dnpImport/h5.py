from .. import dnpdata, create_workspace
import warnings
import numpy as np
import h5py


def saveh5(dataDict, path, overwrite = False):
    '''Save workspace in .h5 format

    Args:
        dataDict (dnpdata_collection): dnpdata_collection object to save.
        path (str): Path to save data
        overwrite (bool): If True, h5 file can be overwritten. Otherwise, h5 file cannot be overwritten
    '''

    if overwrite:
        mode = 'w'
    else:
        mode = 'w-'

    keysList = dataDict.keys()

    f = h5py.File(path, mode)

    for key in keysList:
        dnpDataObject = dataDict[key]

        dnpDataGroup = f.create_group(key, track_order = True)
        if isinstance(dnpDataObject, dnpdata):
            write_dnpdata(dnpDataGroup, dnpDataObject)
        elif isinstance(dnpDataObject, dict):
            write_dict(dnpDataGroup, dnpDataObject)
        else:
            warnings.warn('Could not write key: %s'%str(key))
        
    f.close()

def write_dnpdata(dnpDataGroup, dnpDataObject):
    '''Takes file/group and writes dnpData object to it

    Args:
        dnpDataGroup: h5 group to save data to
        dnpDataObject: dnpdata object to save in h5 format
    '''
    dnpDataGroup.attrs['dnplab_version'] = dnpDataObject.version
    dnpDataGroup.attrs['dnplab_data_type'] = 'dnpdata'
    dims_group = dnpDataGroup.create_group('dims') # dimension names e.g. x,y,z
    attrs_group = dnpDataGroup.create_group('attrs') # dictionary information
    dnp_dataset = dnpDataGroup.create_dataset('values', data = dnpDataObject.values)

    # Save axes information
    for ix in range(len(dnpDataObject.coords)):
        label = dnpDataObject.dims[ix]
        this_axes = dnpDataObject.coords[ix]
        dims_group.create_dataset(label,data = this_axes)
        dims_group[label].make_scale(label)

        dnp_dataset.dims[ix].attach_scale(dims_group[label])

    # Save Parameters
    for key in dnpDataObject.attrs:
        attrs_group.attrs[key] = dnpDataObject.attrs[key]

    # Save proc_steps
    if hasattr(dnpDataObject, 'proc_attrs'):
        proc_attrs = dnpDataObject.proc_attrs
        proc_attrs_group = dnpDataGroup.create_group('proc_attrs', track_order = True)
        for ix in range(len(proc_attrs)):
            proc_step_name = proc_attrs[ix][0]
            proc_dict = proc_attrs[ix][1]
            proc_attrs_group_subgroup = proc_attrs_group.create_group('%i:%s'%(ix,proc_step_name))
            for key in proc_dict:
                value = proc_dict[key]
                proc_attrs_group_subgroup.attrs[key] = value

def write_dict(dnpDataGroup, dnpDataObject):
    '''Writes dictionary to h5 file
    '''
#    dnpDataGroup.attrs['dnplab_version'] = dnpDataObject.version
    dnpDataGroup.attrs['dnplab_data_type'] = 'dict'
#    dnpDataGroup.attrs['dnplab_version'] = dnpDataObject.version
    attrs_group = dnpDataGroup.create_group('attrs')

    for key in dnpDataObject.keys():
        attrs_group.attrs[key] = dnpDataObject[key]


def loadh5(path):
    '''Returns Dictionary of dnpDataObjects

    Args:
        path (str): Path to h5 file

    Returns:
        dnpdata_collection: workspace object with data
    '''

    dnp_dict = {}

    f = h5py.File(path,'r')
    keys_list = f.keys()

    for key in keys_list:

        if f[key].attrs['dnplab_data_type'] == 'dnpdata':
            data = read_dnpdata(f[key])
        elif f[key].attrs['dnplab_data_type'] == 'dict':
            data = read_dict(f[key])
        else:
            warnings.warn('could not import key: %s'%str(key))

        dnp_dict[key] = data
    return create_workspace(dnp_dict)

def read_dnpdata(dnpdata_group):
    axes = []
    dims = []
    attrs = {}
    values = dnpdata_group['values'][:]
    version = dnpdata_group.attrs['dnplab_version']

    for index in range(len(np.shape(values))):
        dim_key = dnpdata_group['values'].dims[index].keys()[0] # assumes 1 key only
        axes.append(dnpdata_group['values'].dims[index][dim_key][:])
        dims.append(dim_key)

    for k in dnpdata_group['attrs'].attrs.keys():
        attrs[k] = dnpdata_group['attrs'].attrs[k]

    data = dnpdata(values, axes, dims, attrs)

    if 'proc_attrs' in dnpdata_group.keys():
        proc_attrs = []
        for k in dnpdata_group['proc_attrs'].keys():
            proc_attrs_name = k.split(':', 1)[1]
            proc_attrs_dict = dict(dnpdata_group['proc_attrs'][k].attrs)
            data.add_proc_attrs(proc_attrs_name, proc_attrs_dict)
    return data

def read_dict(dnpdata_group):
    data = dict(dnpdata_group['attrs'].attrs)
    return data

if __name__ == '__main__':
    pass
