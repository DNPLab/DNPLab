from xml.etree import ElementTree as _ET
import pathlib as _pathlib
import warnings as _warnings
import struct as _struct
import numpy as _np
from .. import DNPData


def import_rs2d(path, datafile="data.dat", headerfile="header.xml", *args, **kwargs):

    path = _pathlib.Path(path)

    #
    # either accept header.xml or data.dat, nothing else for now
    #
    if path.suffix == ".dat":
        path = path.with_name(headerfile)
    if path.suffix != ".xml":
        _warnings.warn(
            "import_rs2d: got file that does not end in .xml, will try to open {} and try to get data.dat".format(
                str(path)
            )
        )

    attrs = _load_rs2d_header(str(path))

    path = path.with_name(datafile)
    data, dims, coords = _load_rs2d_data(path, attrs, **kwargs)

    data = DNPData(data, dims, coords, attrs=attrs)
    dims.reverse()
    data.reorder(dims)
    data.squeeze()

    return data


def _load_rs2d_header(path):

    tree = _ET.parse(path)
    root = tree.getroot()

    attrs = {}
    # params child is important, but for now use all
    for index, child in enumerate(root):
        temp_attrs = {}
        attrs["tag_%i" % (index + 1)] = child.tag
        for ind, entry in enumerate(child):
            try:
                key = entry.find("key").text
                value = entry.find("value").find("value").text
                try:
                    if value.isdigit():
                        value = int(value)
                    else:
                        value = float(value)
                except:
                    continue
                temp_attrs.__setitem__(key, value)
            except AttributeError as e:
                _warnings.warn(
                    "Error in finding key or value at entry {}, skipping entry without elements, attribute error: {}".format(
                        ind, e
                    )
                )

        attrs = {**attrs, **temp_attrs}

    return attrs


def _load_rs2d_data(path, attrs, **kwargs):
    #
    # currently redas whole file in one chunk, you better have enough ram
    #
    endianess = kwargs.get("endianess", ">")
    fmt = kwargs.get("fmt", "f")  # 32bit float
    bitsize = kwargs.get("fmt_size", 4)

    with open(str(path), "rb") as f:
        raw = f.read()

    # change later to better
    assert len(raw) % 2 == 0
    size = int(len(raw) / bitsize)
    data = _np.array(_struct.unpack(endianess + str(size) + fmt, raw))

    data_real = data[slice(0, None, 2)]
    data_imag = data[slice(1, None, 2)]
    data = _np.array(data_real - 1j * data_imag, dtype=complex)
    data *= 1j
    dimNames = list(
        reversed(
            ["ACQUISITION_MATRIX_DIMENSION_" + str(k) + "D" for k in range(1, 5)]
            + ["RECEIVER_COUNT"]
        )
    )
    dims = list(reversed(["t" + str(k) for k in range(len(dimNames))]))
    dimValues = [int(attrs.get(k, 1)) for k in dimNames]

    data = _np.reshape(data, dimValues)

    coords = [_np.arange(k) for k in dimValues]
    try:
        coords[-1] = coords[-1] * float(attrs.get("DWELL_TIME", 1))
        coords[-2] = coords[-2] * float(attrs.get("Polarisation_Growth_Delay", 1))

    except:
        pass

    return data, dims, coords
