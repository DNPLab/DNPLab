import os
from . import *

from ..core.util import concat


def load(path, data_type=None, dim=None, coord=[], verbose=False, *args, **kwargs):
    """Import data from different spectrometer formats

    Args:
        path (str, list): Path to data directory or list of directories
        data_type (str): Type of spectrometer data to import (optional). Allowed values: "prospa", "topspin", "delta", "vnmrj", "tnmr", "specman", "xenon", "xepr", "winepr", "esp", "h5", "power", "vna", "cnsi_powers"
        dim (str): If giving directories as list, name of dimension to concatenate data along
        coord (numpy.ndarray): If giving directories as list, coordinates of new dimension
        verbose (bool): If true, print debugging output
        args: Args passed to spectrometer specific import function
        kwargs: Key word args passed to spectrometer specific import function

    Returns:
        data (dnpData): Data object

    Examples:

        Load a data file

        >>> data = dnp.load('Path/To/File')

        Load a list of files and concatenate down a new dimension called 't1' with coordinates

        >>> data = dnp.load(['1/data.1d','2/data.1d','3/data.1d'], dim = 't1', coord = np.r_[0.1,0.2,0.3])
    """

    if isinstance(path, list):
        if len(coord) != len(path):
            raise ValueError(
                "coord must be a list or array equal in len to the number of paths given"
            )

        data_list = []
        if dim is None:
            dim = "unnamed"
        for filename in path:
            data = load_file(
                filename, data_type=data_type, verbose=verbose, *args, **kwargs
            )
            data_list.append(data)
        # coord could be empty list
        if len(coord) == 0:
            coord = None  # to not break concat call signature

        return concat(data_list, dim=dim, coord=coord)

    else:
        return load_file(path, data_type=data_type, verbose=verbose, *args, **kwargs)


def load_file(path, data_type=None, verbose=False, *args, **kwargs):
    """Import data from different spectrometer formats

    Args:
        path (str): Path to data directory or file
        data_type (str): Type of spectrometer data to import (optional). Allowed values: "prospa", "topspin", "delta", "vnmrj", "tnmr", "specman", "xenon", "xepr", "winepr", "esp", "h5", "power", "vna", "cnsi_powers"
        verbose (bool): If true, print additional debug outputs
        args: Arguments passed to spectrometer specific import function
        kwargs: Key word arguments passed to spectrometer specific import function

    Returns:
        data (dnpData): Data object
    """

    path = os.path.normpath(path)
    if os.path.isdir(path) and path[-1] != os.sep:
        path = path + os.sep

    if data_type == None:
        data_type = autodetect(path, verbose=verbose)

    if data_type == "prospa":
        return prospa.import_prospa(path, *args, **kwargs)

    elif data_type == "topspin":
        return topspin.import_topspin(path, verbose=verbose, *args, **kwargs)

    elif data_type == "topspin pdata":
        # import_topspin should also handle this type, this is a workaround
        return topspin.load_pdata(path, verbose=verbose, *args, **kwargs)

    elif data_type == "topspin dir":
        return topspin.import_topspin_dir(path, *args, **kwargs)

    elif data_type == "delta":
        return delta.import_delta(path, *args, **kwargs)

    elif data_type == "vnmrj":
        return vnmrj.import_vnmrj(path, *args, **kwargs)

    elif data_type == "tnmr":
        return tnmr.import_tnmr(path, *args, **kwargs)

    elif data_type == "specman":
        return specman.import_specman(path, *args, **kwargs)

    elif data_type in ["xepr", "xenon"]:
        return bes3t.import_bes3t(path, *args, **kwargs)

    elif data_type in ["winepr", "esp"]:
        return winepr.import_winepr(path, *args, **kwargs)

    elif data_type == "h5":
        return h5.load_h5(path, *args, **kwargs)

    elif data_type == "power":
        return power.importPower(path, *args, **kwargs)

    elif data_type == "vna":
        return vna.import_vna(path, *args, **kwargs)

    elif data_type == "cnsi_powers":
        return cnsi.get_powers(path, *args, **kwargs)

    else:
        raise ValueError("Invalid data type: %s" % data_type)


# TODO rename to detect_file_format
def autodetect(test_path, verbose=False):
    """Automatically detect spectrometer format

    Args:
        test_path (str): Test directory
        verbose (bool): If true, print output for debugging

    Returns:
        str: Spectrometer type as string

    """

    if verbose == True:
        print("current directory:", os.getcwd())
        print("data path:", test_path)
        abs_path = os.path.abspath(test_path)
        print("absolute path:", abs_path)

    # Remove trailing separator
    if test_path[-1] == os.sep:
        test_path = test_path[:-1]
        if verbose:
            print("removed trailing separator:", os.sep)

    path_exten = os.path.splitext(test_path)[1]
    if path_exten != "" and verbose:
        print("Extension:", path_exten)

    if path_exten == ".DSC" or path_exten == ".DTA" or path_exten == ".YGF":
        spectrometer_format = "xepr"
    elif path_exten in [".par", ".spc"]:
        spectrometer_format = "winepr"
    elif path_exten in [".d01", ".exp"]:
        spectrometer_format = "specman"
    elif path_exten == ".jdf":
        spectrometer_format = "delta"
    elif (
        os.path.isdir(test_path)
        #        and ("fid" in os.listdir(test_path) or "ser" in os.listdir(test_path))
        and ("acqu" in os.listdir(test_path) or "acqus" in os.listdir(test_path))
    ):
        spectrometer_format = "topspin"
    elif os.path.isdir(test_path) and (
        "proc" in os.listdir(test_path) or "procss" in os.listdir(test_path)
    ):
        spectrometer_format = "topspin pdata"
    elif os.path.isdir(test_path) and path_exten == ".fid":
        spectrometer_format = "vnmrj"
    elif path_exten in [".1d", ".2d", ".3d", ".4d"]:
        spectrometer_format = "prospa"
    elif (
        os.path.isdir(test_path)
        and "acqu.par" in os.listdir(test_path)
        and "data.csv" in os.listdir(test_path)
    ):
        spectrometer_format = "prospa"
    elif path_exten == ".h5":
        spectrometer_format = "h5"
    else:
        raise TypeError(
            "No data type given and autodetect failed to detect format, please specify a format"
        )

    if verbose:
        print("Spectrometer Format:", spectrometer_format)

    return spectrometer_format
