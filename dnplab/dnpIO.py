from . import dnpImport


def import_data(path, data_type=None, *args, **kwargs):
    """Import data from different spectrometer formats

    Args:
        path (str): Path to data directory or file
        data_type (str): Type of data to import

    Returns:
        data (dnpData): Data object
    """

    if data_type == None:
        data_type = autodetect(path)

    if data_type == NotImplemented:
        return NotImplemented

    if data_type == "prospa":
        return dnpImport.prospa.import_prospa(path, *args, **kwargs)

    elif data_type == "topspin":
        return dnpImport.topspin.import_topspin(path, *args, **kwargs)

    elif data_type == "topspin dir":
        return dnpImport.topspin.import_topspin_dir(path, *args, **kwargs)

    elif data_type == "vnmrj":
        return dnpImport.vnmrj.import_vnmrj(path, *args, **kwargs)


def autodetect(path):
    """Automatically detect type of data in directory"""

    return NotImplemented
