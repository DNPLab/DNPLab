import sys
import numpy as _np
sys.path.insert(0,r"C:\Users\krieger\Personal\MRD_Projects\DNPLab")
import dnplab as dnp

def signal_to_noise(
    data: dnp.DNPData,
    signal_region: list = slice(0, None),
    noise_region: list = (None, None),
    dim: str = "f2",
    remove_background: list = None,
    **kwargs
):
    """Find signal-to-noise ratio

    Simplest implementation: select largest value in a signal_region and divide this value by the estimated std. deviation of another noise_region. If the noise_region list contains (None,None) (the default) then all points except the points +10% and -10% around the maximum are used for the noise_region.

    Args:
        data: Spectrum data
        signal_region (list): list with a single tuple (start,stop) of a region where a signal should be searched, default is [slice(0,None)] which is the whole spectrum
        noise_region (list): list with tuples (start,stop) of regions that should be taken as noise, default is (None,None)
        dim (str): dimension of data that is used for snr calculation, default is 'f2'
        remove_background (list): if this is not None (a list of tuples, or a single tuple) this will be forwarded to dnp.remove_background, together with any kwargs
        kwargs : parameters for dnp.remove_background

    Returns:
        SNR (float): Signal to noise ratio

    Examples:

        A note for the usage: regions can be provided as (min,max), slices use indices.
        To use the standard values just use

            >>> snr = dnp.signal_to_noise(data)

        If you want to select a region for the noise and the signal:

            >>> snr = dnp.signal_to_noise(data,[(-1.23,300.4)],noise_region=[(-400,-240.5),(123.4,213.5)])

        With background subtracted:

            >>> snr = dnp.signal_to_noise(data,[(-1.23,300.4)],noise_region=[(-400,-240.5),(123.4,213.5)],remove_background=[(123.4,213.5)])

        This function allows to use a single tuple instead of a list with a single tuple for signal_region, noise_region and remove_background. This is for convenience, slices are currently only supoprted for signal_region and noise_region.

            >>> snr = dnp.signal_to_noise(data,(-1.23,300.4),noise_region=[(-400,-240.5),(123.4,213.5],remove_background=(123.4,213.5))

    """
    import warnings
    import scipy.optimize as _scipy_optimize

    # convenience for signal and noise region
    def _convenience_tuple_to_list(possible_region: list):
        if possible_region is None:
            return possible_region
        # we assume its iterable
        try:
            l = len(possible_region)
            if l != 2:
                return possible_region  # it seems to be iterable?
        except TypeError:
            return [possible_region]  # return as is in a list, might be slice
        try:
            # check whether we can interpret it as value
            a = int(possible_region[0])
            return [
                (possible_region[0], possible_region[1])
            ]  # make a list that contains a tuple
        except TypeError:
            if type(possible_region) == list:
                return possible_region
            return [possible_region]

    if not (dim in data.dims):
        raise ValueError("dim {0} not in data.dims ({1})".format(dim, data.dims))

    #multiple dimensions behaviour:
    if len(data.shape)==2:
        ind=data.index(dim)



    # remove background
    remove_background = _convenience_tuple_to_list(remove_background)
    if remove_background is not None:
        deg = kwargs.pop("deg", 1)
        data = dnp_remove_background(data, dim, deg, remove_background)

    # concat noise region to one region
    if type(noise_region)==float:
        noise = noise_region
    else:
        noise_region = _convenience_tuple_to_list(noise_region)
        if (None, None) in noise_region:
            signal_arg = _np.argmax(_np.abs(data[dim, signal_region[0]]))
            datasize = data[dim, :].size
            noise_region = [
                slice(0, int(_np.maximum(2, int(signal_arg * 0.9)))),
                slice(int(_np.minimum(datasize - 2, int(signal_arg * 1.1))), None),
            ]
        # concatenate noise_regions
        noise_0 = _np.abs(data[dim, noise_region[0]])

        for k in noise_region[1:]:
            noise_0.concatenate(_np.abs(data[dim, k]), dim)

        noise = float( _np.std(_np.abs(noise_0[dim, slice(0, None)])) )

    # currently only one method avaiable -> absolute value
    signal_region = _convenience_tuple_to_list(signal_region)
    if len(signal_region) > 1:

        snr = [_np.max(_np.abs(data[dim, signal_region[0]])) ]
        added_snr = signal_to_noise(data,signal_region=signal_region[1:],noise_region=noise)
        if type(added_snr) == list:
            snr += added_snr
        else:
            snr += [added_snr]
        return np.array(snr)

    signal = _np.max(_np.abs(data[dim, signal_region[0]]))

    return float(signal / noise)


if __name__=="__main__":

    # snr dnplab implementation
    #example data:
    data = dnp.load("../../data/prospa/water_phase_cycled/data.2d")

    # data has now f2 and Average as dimensions
    data = dnp.fourier_transform(data)

    coords=data.coords['f2']
    argmax=_np.argmax(data,axis='f2')
    print(data.shape)

    """
    SNR change, possible list of snr_regions
    - 1d -> take axis
    - 2d -> make automatically 2d snr values
    - 3d or more: define possible second dimension and if not raise error
    - accept list of snr regions
     - accept list of at least 1 noise region
    """
    print(data.dims)

    print(argmax,coords[0:10])
