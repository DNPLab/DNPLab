import numpy as np

def remove_offset(data, dim="t2", offset_points=10):
    """Remove DC offset from FID by averaging the last few data points and subtracting the average

    Args:
        all_data (dnpdata, dict): Data container for data

    +---------------+------+---------+----------------------------------------------------------+
    | parameter     | type | default | description                                              |
    +===============+======+=========+==========================================================+
    | dim           | str  | 't2'    | Dimension to calculate DC offset                         |
    +---------------+------+---------+----------------------------------------------------------+
    | offset_points | int  | 10      | Number of points at end of data to average for DC offset |
    +---------------+------+---------+----------------------------------------------------------+

    Returns:
        dnpdata: data object with offset removed
    """

    proc_parameters = {
        "dim": dim,
        "offset_points": offset_points,
    }

    dim = proc_parameters["dim"]
    offset_points = int(proc_parameters["offset_points"])

    offsetData = data[dim, -1 * offset_points :].values
    offsetData = offsetData.reshape(-1)
    offset = np.mean(offsetData)

    data -= offset

    proc_attr_name = "remove_offset"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    return data


def background(data, dim = "t2", deg = 0, regions = None):
    """ Remove background from data
    """

    out = data.copy()
    out.unfold(dim)

    coord = out.coords[dim]
    fit_points = [False for x in coord]

    for region in regions:
        fit_points = [fit_points[ix] or ((coord[ix] >= region[0]) & (coord[ix] <= region[1])) for ix in range(len(coord))]


    for ix in range(out.shape[1]):
        p = np.polyfit(coord[fit_points], out.values[:,ix][fit_points], deg = deg)
        fit = np.polyval(p,coord)
        out.values[:,ix] = fit

    out.fold()

    return out


