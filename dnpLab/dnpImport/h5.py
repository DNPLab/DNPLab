from .. import dnpData
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
        dnpDataGroup.attrs['dnpLab_version'] = dataDict[key].version
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
        if 'proc_steps' in dnpDataObject.__dict__:
            proc_steps = dnpDataObject.proc_steps
            proc_steps_group = dnpDataGroup.create_group('proc_steps', track_order = True)
            for ix in range(len(proc_steps)):
                proc_step_name = proc_steps[ix][0]
                proc_dict = proc_steps[ix][1]
                proc_steps_group_subgroup = proc_steps_group.create_group('%i:%s'%(ix,proc_step_name))
                for key in proc_dict:
                    value = proc_dict[key]
                    proc_steps_group_subgroup.attrs[key] = value
#            proc_steps_group.attrs['proc_steps'] = dnpDataObject.proc_steps

    f.close()

def loadh5(path):
    '''
    Returns Dictionary of dnpDataObjects
    '''

    dnpDict = {}

    f = h5py.File(path,'r')
    keysList = f.keys()
#    print('keys:')
#    print(keysList)

    for key in keysList:
        axes = []
        dims = []
        attrs = {}
        data = f[key]['values'][:]
        version = f[key].attrs['dnpLab_version']

        for index in range(len(np.shape(data))):
            dimKey = f[key]['values'].dims[index].keys()[0] # assumes 1 key only
            axes.append(f[key]['values'].dims[index][dimKey][:])
            dims.append(dimKey)

        for k in f[key]['attrs'].attrs.keys():
#            print(k)
#            print(f[key]['attrs'].attrs[k])
            attrs[k] = f[key]['attrs'].attrs[k]

        dnpDict[key] = dnpData(data,axes,dims,attrs)

        if 'proc_steps' in f[key].keys():
            proc_steps = []
            for k in f[key]['proc_steps'].keys():
                proc_step_name = k.split(':', 1)[1]
                proc_step_dict = dict(f[key]['proc_steps'][k].attrs)
                dnpDict[key].add_proc_step(proc_step_name, proc_step_dict)



    return dnpDict

