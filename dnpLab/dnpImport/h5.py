from .. import dnpData, create_workspace
import numpy as np
import h5py


def saveh5(dataDict, path, overwrite = False):
    '''
    Save All Data in .h5 format
    Requires a dictionary of dnpData objects

    '''

    if overwrite:
        mode = 'w'
    else:
        mode = 'w-'

    keysList = dataDict.keys()

    f = h5py.File(path, mode)

#    f.attrs['dnpLab_version'] = version

    for key in keysList:
        dnpDataObject = dataDict[key]
        dnpDataGroup = f.create_group(key, track_order = True)

        write_dnpdata(dnpDataGroup, dnpDataObject)
        
    f.close()

def write_dnpdata(dnpDataGroup, dnpDataObject):
    '''Takes file/group and writes dnpData object to it
    '''
#    dnpDataGroup.attrs['dnpLab_version'] = dataDict[key].version
    dnpDataGroup.attrs['dnpLab_version'] = dnpDataObject.version
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
    if 'proc_attrs' in dnpDataObject.__dict__:
        proc_attrs = dnpDataObject.proc_attrs
        proc_attrs_group = dnpDataGroup.create_group('proc_attrs', track_order = True)
        for ix in range(len(proc_attrs)):
            proc_step_name = proc_attrs[ix][0]
            proc_dict = proc_attrs[ix][1]
            proc_attrs_group_subgroup = proc_attrs_group.create_group('%i:%s'%(ix,proc_step_name))
            for key in proc_dict:
                value = proc_dict[key]
                proc_attrs_group_subgroup.attrs[key] = value

def loadh5(path):
    '''
    Returns Dictionary of dnpDataObjects
    '''

    dnp_dict = {}

    f = h5py.File(path,'r')
    keys_list = f.keys()
#    print('keys:')
#    print(keysList)

    for key in keys_list:
        axes = []
        dims = []
        attrs = {}
        data = f[key]['values'][:]
        version = f[key].attrs['dnpLab_version']

        for index in range(len(np.shape(data))):
            dim_key = f[key]['values'].dims[index].keys()[0] # assumes 1 key only
            axes.append(f[key]['values'].dims[index][dim_key][:])
            dims.append(dim_key)

        for k in f[key]['attrs'].attrs.keys():
#            print(k)
#            print(f[key]['attrs'].attrs[k])
            attrs[k] = f[key]['attrs'].attrs[k]

        dnp_dict[key] = dnpData(data,axes,dims,attrs)

        if 'proc_attrs' in f[key].keys():
            proc_attrs = []
            for k in f[key]['proc_attrs'].keys():
                proc_attrs_name = k.split(':', 1)[1]
                proc_attrs_dict = dict(f[key]['proc_attrs'][k].attrs)
                dnp_dict[key].add_proc_attrs(proc_attrs_name, proc_attrs_dict)



    return create_workspace(dnp_dict)

