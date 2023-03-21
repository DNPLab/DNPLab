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

# config for load_csv file, can be adjusted to add more options (dialect?)
# could also be checked for userconfig, bt this is not implemented right now
parent_dir = pathlib.Path(pathlib.Path(__file__).absolute().parents[1])
p0 = parent_dir.joinpath("config/io_config.conf")
config_io = configparser.ConfigParser()
config_io.read(str(p0))


def load_csv(
    filename,
    tcol=0,
    real=1,
    imag=2,
    skiprows=0,
    maxrows=-1,
    convert_time=lambda x: float(x.replace(",", ".")) / 1e6,
    convert_data=lambda x: float(x.replace(",", ".")),
    **kwargs
):
    """
    load_csv files function

    Args:
        filename: str/path like
        tcol: column index for time data (default=0), None = not applicable, will count from 0 to npoints and then apply covert_time!
        real: column index for real part (default=1), None = not applicable, will be set to zero
        imag: column index for imaginary part (default=2), None = not applicable, will be set to zero
        skiprows: number of rows to skip at beginning (default=0)
        maxrows: if this is larger than zeros read at most maxrows rows
        convert_time: callable that converts the time strings to a number (default: assumes us and converts into s)
        convert_data: callable that converts data to a number (default: replaces , with .)

        delimiter: optional, sets the delimiter in the csv file (default ; ), can be set in ../config/io_config.conf
        dims: optional, sets name for dimension (default: t2)
        **kwargs are forwarded to csv.reader object


    Returns:
        DNPData: data object with values,coords and dim


    example:

        data=load_csv('csv_arrLNA_data.csv',imag=2,skiprows=1)
    """

    # standards from config (?)
    delimiter = config_io.get(
        "load_csv", "delimiter"
    )  # or config_io['load_csv']['delimiter']

    if "delimiter" in kwargs.keys():
        delimiter = kwargs.pop("delimiter")

    dims = kwargs.pop("dims", "t2")

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
