import os
from . import dnpIO
import os


def load(path, data_type=None, *args, **kwargs):
    """Import data from different spectrometer formats
    +-----------------+
    | Supported Types |
    +-----------------+
    | prospa          |
    +-----------------+
    | topspin         |
    +-----------------+
    | vnmrj           |
    +-----------------+
    | specman         |
    +-----------------+
    Args:
        path (str): Path to data directory or file
        data_type (str): Type of spectrometer data to import
    Returns:
        data (dnpData): Data object
    """
    path = os.path.normpath(path)
    if os.path.isdir(path) and path[-1] != os.sep:
        path = path + os.sep

    if data_type == None:
        data_type = autodetect(path)

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

    elif data_type == "h5":
        return dnpIO.h5.loadh5(path, *args, **kwargs)

    elif data_type == "power":
        return dnpIO.power.importPower(path, *args, **kwargs)

    elif data_type == "vna":
        return dnpIO.vna.import_vna(path, *args, **kwargs)

    else:
        raise ValueError("Invalid data type: %s" % data_type)


def autodetect(path):
    """Automatically detect type of data in directory"""

    is_file = os.path.isfile(path)

    # If path is file
    if is_file:
        path, filename = os.path.split(path)
        _, extension = os.path.splitext(filename)

        if extension == '.h5':
            return 'h5'
        elif extension in ['.1d', '.2d', '.3d', '.4d']:
            return 'prospa'
        elif filename == 'fid':
            files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
            if 'procpar' in files:
                return 'vnmrj'
            else:
                return 'topspin'

    # if path is directory
    else:
        if path[-4:] == '.fid':
            return 'vnmrj'

        files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]

        if ('fid' in files) or ('ser' in files):
            return 'topspin'

        elif 'acqu.par' in files:
            return 'prospa'

