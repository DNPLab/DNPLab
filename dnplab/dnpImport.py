import os
from . import dnpIO


def load(path, data_type=None, *args, **kwargs):
    """Import data from different spectrometer formats
    +-----------------+
    | Supported Types |
    +-----------------+
    | prospa          |
    +-----------------+
    | topspin         |
    +-----------------+
    | delta           |
    +-----------------+
    | vnmrj           |
    +-----------------+
    | specman         |
    +-----------------+
    | xenon and xepr  |
    +-----------------+
    | winepr and esp  |
    +-----------------+
    | h5              |
    +-----------------+
    | power           |
    +-----------------+
    | vna             |
    +-----------------+
    | cnsi_powers     |
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

    elif data_type == "delta":
        return dnpIO.delta.import_delta(path, *args, **kwargs)

    elif data_type == "vnmrj":
        return dnpIO.vnmrj.import_vnmrj(path, *args, **kwargs)

    elif data_type == "specman":
        return dnpIO.specman.import_specman(path, *args, **kwargs)

    elif data_type == "xepr" or data_type == "xenon":
        return dnpIO.bes3t.import_bes3t(path, *args, **kwargs)

    elif data_type == "winepr" or data_type == "esp":
        return dnpIO.parspc.import_parspc(path, *args, **kwargs)

    elif data_type == "h5":
        return dnpIO.loadh5.load_h5(path, *args, **kwargs)

    elif data_type == "power":
        return dnpIO.power.importPower(path, *args, **kwargs)

    elif data_type == "vna":
        return dnpIO.vna.import_vna(path, *args, **kwargs)

    elif data_type == "cnsi_powers":
        return dnpIO.cnsi.get_powers(path, *args, **kwargs)

    else:
        raise ValueError("Invalid data type: %s" % data_type)


def autodetect(test_path):
    """Automatically detect type of data in directory"""

    if test_path[-1] == os.sep:
        test_path = test_path[:-1]

    if test_path[-4:] == ".DSC" or test_path[-4:] == ".DTA" or test_path[-4:] == ".YGF":
        type = "xepr"
    elif test_path[-4:] == ".par" or test_path[-4:] == ".spc":
        type = "winepr"
    elif test_path[-4:] == ".d01" or test_path[-4:] == ".exp":
        type = "specman"
    elif test_path[-4:] == ".jdf":
        type = "delta"
    elif os.path.isdir(test_path) and "pdata" in os.listdir(test_path):
        type = "topspin"
    elif os.path.isdir(test_path) and test_path[-4:] == ".fid":
        type = "vnmrj"
    elif (
        os.path.isdir(test_path)
        and "acqu.par" in os.listdir(test_path)
        and "data.csv" in os.listdir(test_path)
    ):
        type = "prospa"
    elif test_path[-3:] == ".h5":
        type = "h5"
    else:
        raise TypeError(
            "No data type given and autodetect failed to detect format, please specify a format"
        )

    return type
