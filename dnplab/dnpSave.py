import os
from . import io, dnpdata, dnpdata_collection, create_workspace


def save(data_object, filename, save_type=None, *args, **kwargs):
    """Save data to h5 format

    Args:
        data_object (dnpdata object) : dnpdata object to save
        filename (str): name of file, must include extension .h5
        save_type (str): Type of file to save (optional). Allowed values: "h5"

    """

    if save_type == None:
        save_type = autodetect(filename)

    if save_type == "h5":
        if isinstance(data_object, dnpdata_collection):
            return dnpIO.h5.save_h5(data_object, filename, *args, **kwargs)
        elif isinstance(data_object, dnpdata):
            ws = create_workspace()
            ws.add("data", data_object)
            return dnpIO.h5.save_h5(ws, filename, *args, **kwargs)
        else:
            raise TypeError(
                "object format not recognized, must be dnpdata or dnpdata_collection"
            )

    else:
        raise TypeError("File type not recognized, you must specify a save format")


def autodetect(test_name):

    if test_name[-1] == os.sep:
        test_name = test_name[:-1]

    if os.path.splitext(test_name)[1] == ".h5":
        type = "h5"
    else:
        raise TypeError("File type not recognized, you must specify a save format")

    return type
