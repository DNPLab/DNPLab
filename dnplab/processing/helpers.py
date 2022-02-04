import numpy as np
from scipy.signal import savgol_filter

from ..core.data import DNPData
from ..processing.integration import integrate


def calculate_enhancement():
    """Calculate enhancement 
    
    Returns:
        NotImplemented

    """
    

    return NotImplemented


def signal_to_noise():
    """Find signal-to-noise ratio

    Returns:
        NotImplemented
    """

    return NotImplemented


def smooth(data, dim="t2", window_length=11, polyorder=3):
    out = data.copy()

    out.unfold(dim)

    out.values = savgol_filter(out.values, window_length, polyorder, axis=0)

    out.fold()

    return out


def left_shift(data, dim="t2", shift_points=0):
    """Remove points from the left of data

    Args:
        all_data (dnpdata, dict): Data container for data

    +---------------+------+---------+--------------------------------------------------+
    | parameter     | type | default | description                                      |
    +===============+======+=========+==================================================+
    | dim           | str  | "t2"    | dimension to shift                               |
    +---------------+------+---------+--------------------------------------------------+
    | shift_points  | int  | 0       | Number of points to remove from left of data     |
    +---------------+------+---------+--------------------------------------------------+

    Returns:
        dnpdata: data object with left-shifted data
    """

    data = data[dim, shift_points:]

    proc_attr_name = "left_shift"
    proc_parameters = {
        "dim": dim,
        "points": shift_points,
    }
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    return data


def normalize():
    return NotImplemented


def reference():
    return NotImplemented
