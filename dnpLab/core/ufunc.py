import defaults
import nddata
import numpy as np

def fourier_transform(data, proc_parameters):
    '''Perform Fourier Transform down dim dimension given in proc_parameters

    .. Note::
        Assumes dt = t[1] - t[0]

    Args:
        data (nddata): Data container
        proc_parameters (dict, procParam): Processing parameters

    Returns:
        nddata: Fourier Transformed data
        
    Example:

    .. code-block:: python
    
        proc_parameters['dim'] = 't'
        proc_parameters['zero_fill_factor'] = 2
        proc_parameters['shift'] = True
        proc_parameters['convert_to_ppm'] = True

        all_data = dnplab.dnpNMR.fourier_transform(all_data, proc_parameters)
    '''

    required_parameters = defaults._fourier_transform
#    proc_parameters = update_parameters(proc_parameters,requiredList,_default_fourier_transform_parameters)

    # Add required parameters to proc_parameters
    print(required_parameters)
    for key in required_parameters:
        if key not in proc_parameters:
            proc_parameters[key] = required_parameters[key]
#
    dim = proc_parameters['dim']
    zero_fill_factor = proc_parameters['zero_fill_factor']
    shift = proc_parameters['shift']
    convert_to_ppm = proc_parameters['convert_to_ppm']

    index = data.dims.index(dim)
    dt = data.coords[index][1] - data.coords[index][0]

    n_pts = zero_fill_factor*len(data.coords[index])
    f = (1./(n_pts*dt))*np.r_[0:n_pts]

    if shift == True:
        f -= (1./(2*dt))

#    if convert_to_ppm:
#        nmr_frequency = data.attrs['nmr_frequency']
#        f /= (nmr_frequency / 1.e6)

    data.values = np.fft.fft(data.values,n=n_pts,axis=index)
    if shift:
        data.values = np.fft.fftshift(data.values,axes=index)
    data.coords[index] = f

    return data


if __name__ == '__main__':
    x = np.r_[0:10]
    y = np.r_[0:20]
    z = np.r_[0:15]
    data = nddata.nddata_core(np.array(range(len(x)*len(y)*len(z))).reshape(len(x),len(y),len(z)), ['x','y','z'], [x, y, z])

    out = fourier_transform(data,{'dim':'x'})
    print(out)

