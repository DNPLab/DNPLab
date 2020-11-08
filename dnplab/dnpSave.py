import os
from . import dnpIO


def save(data_object, filename, save_type=None, *args, **kwargs):
    """Save data to h5 format
    +-----------------+
    | Supported Types |
    +-----------------+
    | h5              |
    +-----------------+

    Args:
        data_object (dnpdata object) : dnpdata object to save
        filename (str): name of saved .csv or .h5 file
        save_type (str): Type of spectrometer data to import

    """

    if save_type == None:
        save_type = autodetect(filename)

    if save_type == "h5":
        return dnpIO.saveh5.save_h5(data_object, filename, *args, **kwargs)

    else:
        raise TypeError("File type not recognized, you must specify a save format")


def autodetect(test_name):
    """Automatically detect type of file to save"""

    if test_name[-1] == os.sep:
        test_name = test_name[:-1]

    if test_name[-3:] == ".h5":
        type = "h5"
    else:
        raise TypeError("File type not recognized, you must specify a save format")

    return type
