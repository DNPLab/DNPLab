import os
from . import dnpIO


def load(path, data_type=None, *args, **kwargs):
    """Import data from different spectrometer formats

    Args:
        path (str): Path to data directory or file
        data_type (str): Type of spectrometer data to import

    Returns:
        data (dnpData): Data object
    """
    if path[-1:] != os.sep:
        path = path + os.sep

    if data_type == None:
        data_type = autodetect(path)

    if data_type == NotImplemented:
        return ValueError("Autodetecting data type not implemented")

    if data_type == "prospa":
        return dnpIO.prospa.import_prospa(path, *args, **kwargs)

    elif data_type == "topspin":
        return dnpIO.topspin.import_topspin(path, *args, **kwargs)

    elif data_type == "topspin dir":
        return dnpIO.topspin.import_topspin_dir(path, *args, **kwargs)

    elif data_type == "vnmrj":
        return dnpIO.vnmrj.import_vnmrj(path, *args, **kwargs)

    elif data_type == "specman":
        return dnpIO.specman.import_specman(path, *args, **kwargs)

    else:
        raise ValueError("Invalid data type: %s" % data_type)


def autodetect(path):
    """Automatically detect type of data in directory"""

    return NotImplemented
