import logging
import csv
import pathlib
import sys
import configparser

import numpy as _np

from .. import DNPData

# define logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler(stream=sys.stdout)
handler.setLevel(logging.INFO)
logger.addHandler(handler)


def load_csv(
    filename,
    tcol=0,
    real=1,
    imag=2,
    skiprows=0,
    maxrows=-1,
    convert_time=lambda x: float(x.replace(",", ".")),
    convert_data=lambda x: float(x.replace(",", ".")),
    **kwargs
):
    """function that loads load_csv files

    Args:
        filename (str): String or path like file that is read
        tcol (int): column index for time data (default=0), None = not applicable, will count from 0 to npoints and then apply covert_time!
        real (int): column index for real part (default=1), None = not applicable, will be set to zero
        imag (int): column index for imaginary part (default=2), None = not applicable, will be set to zero
        skiprows (int): number of rows to skip at beginning (default=0)
        maxrows (int): if this is larger than zero, read at most maxrows rows. (default: -1)
        convert_time (callable): callable that converts the time strings to a number (default: replaces , with .)
        convert_data (callable): callable that converts data to a number (default: replaces , with .)

        delimiter (str): optional, sets the delimiter in the csv file (default ; )
        dims (str,list): optional, sets name for dimension (default: ["t2"])
        **kwargs are forwarded to csv.reader object


    Returns:
       data (dnpData): Data object with values,coords and dim


    example:

        data=load_csv('csv_arrLNA_data.csv',imag=2,skiprows=1)
    """

    delimiter = kwargs.pop("delimiter", ";")

    dims = kwargs.pop("dims", ["t2"])

    def _checknone(x, row, ind=None):
        if x is None and (ind is not None):
            return str(ind)
        if x is None:
            return str(0)
        return row[x]

    data = []
    with open(filename, newline="") as csvfile:
        reader = csv.reader(csvfile, delimiter=delimiter, **kwargs)
        for row_ind, row in enumerate(reader):
            if row_ind < skiprows:
                continue
            if row_ind > maxrows and maxrows > 0:
                break
            try:
                t = convert_time(_checknone(tcol, row, row_ind - skiprows))
                r = convert_data(_checknone(real, row))
                i = convert_data(_checknone(imag, row))
            except IndexError:
                raise IndexError(
                    "Index Error while accesing csv rows, one of tcol/real/imag ({0}/{1}/{2}) is out of range or is invalid".format(
                        tcol, real, imag
                    )
                )

            data.append([t, r + 1j * i])

    data = _np.array(data, dtype=complex)
    values = data[:, 1]
    coords = _np.real(data[:, 0])

    logger.debug(
        "Loading data with {0} elements and from file {1} with delimiter {2}".format(
            len(values), filename, delimiter
        )
    )

    return DNPData(values, dims, [coords])
