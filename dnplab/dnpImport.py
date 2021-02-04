import os
from . import dnpIO


def load(path, data_type=None, *args, **kwargs):
    """Import data from different spectrometer formats

    +--------------------+--------------+--------------------+
    | parameter          | type         | allowed values     |
    +====================+==============+====================+
    | data_type          | str          | "prospa"           |
    +--------------------+--------------+--------------------+
    |                    |              | "topspin"          |
    +--------------------+--------------+--------------------+
    |                    |              | "delta"            |
    +--------------------+--------------+--------------------+
    |                    |              | "vnmrj"            |
    +--------------------+--------------+--------------------+
    |                    |              | "tnmr"             |
    +--------------------+--------------+--------------------+
    |                    |              | "specman"          |
    +--------------------+--------------+--------------------+
    |                    |              | "xenon" or "xepr"  |
    +--------------------+--------------+--------------------+
    |                    |              | "winepr" or "esp"  |
    +--------------------+--------------+--------------------+
    |                    |              | "h5"               |
    +--------------------+--------------+--------------------+
    |                    |              | "power"            |
    +--------------------+--------------+--------------------+
    |                    |              | "vna"              |
    +--------------------+--------------+--------------------+
    |                    |              | "cnsi_powers"      |
    +--------------------+--------------+--------------------+

    Args:
        path (str): Path to data directory or file
        data_type (str): Type of spectrometer data to import (optional)

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

    elif data_type == "tnmr":
        return dnpIO.tnmr.import_tnmr(path, *args, **kwargs)

    elif data_type == "specman":
        return dnpIO.specman.import_specman(path, *args, **kwargs)

    elif data_type == "xepr" or data_type == "xenon":
        return dnpIO.bes3t.import_bes3t(path, *args, **kwargs)

    elif data_type == "winepr" or data_type == "esp":
        return dnpIO.winepr.import_winepr(path, *args, **kwargs)

    elif data_type == "h5":
        return dnpIO.h5.load_h5(path, *args, **kwargs)

    elif data_type == "power":
        return dnpIO.power.importPower(path, *args, **kwargs)

    elif data_type == "vna":
        return dnpIO.vna.import_vna(path, *args, **kwargs)

    elif data_type == "cnsi_powers":
        return dnpIO.cnsi.get_powers(path, *args, **kwargs)

    else:
        raise ValueError("Invalid data type: %s" % data_type)


def autodetect(test_path):

    if test_path[-1] == os.sep:
        test_path = test_path[:-1]

    path_exten = os.path.splitext(test_path)[1]
    if path_exten == ".DSC" or path_exten == ".DTA" or path_exten == ".YGF":
        type = "xepr"
    elif path_exten == ".par" or path_exten == ".spc":
        type = "winepr"
    elif path_exten == ".d01" or path_exten == ".exp":
        type = "specman"
    elif path_exten == ".jdf":
        type = "delta"
    elif (
        os.path.isdir(test_path)
        and "pdata" in os.listdir(test_path)
        and "acqus" in os.listdir(test_path)
    ):
        type = "topspin"
    elif os.path.isdir(test_path) and path_exten == ".fid":
        type = "vnmrj"
    elif path_exten in [".1d", ".2d", ".3d", ".4d"]:
        type = "prospa"
    elif (
        os.path.isdir(test_path)
        and "acqu.par" in os.listdir(test_path)
        and "data.csv" in os.listdir(test_path)
    ):
        type = "prospa"
    elif path_exten == ".h5":
        type = "h5"
    else:
        raise TypeError(
            "No data type given and autodetect failed to detect format, please specify a format"
        )

    return type
