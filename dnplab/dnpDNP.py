from . import dnpNMR, dnpdata, dnpdata_collection
import numpy as np
import copy


def return_data(all_data):

    is_workspace = False
    if isinstance(all_data, dnpdata):
        data = all_data.copy()
    elif isinstance(all_data, dict):
        raise ValueError("Type dict is not supported")
    elif isinstance(all_data, dnpdata_collection):
        is_workspace = True
        if all_data.processing_buffer in all_data.keys():
            data = all_data[all_data.processing_buffer]
        else:
            raise ValueError("No data in processing buffer")
    else:
        raise ValueError("Data type not supported")

    return data, is_workspace


def calculate_enhancement(
    all_data,
    off_spectrum=1,
    on_spectra="all",
    integrate_center=0,
    integrate_width="full",
    method="integrate",
    dim="t2",
):
    """Calculate enhancement from DNP data

    Args:
        all_data (dnpdata, dict): data container

    +------------------+----------------------------+-------------+----------------------------------------------------------------------+
    | parameter        | type                       | default     | description                                                          |
    +------------------+----------------------------+-------------+----------------------------------------------------------------------+
    | off_spectrum     | int or dnpdata             | 1           | slice of 2D data to be used as p = 0 spectrum                        |
    +------------------+----------------------------+-------------+----------------------------------------------------------------------+
    | on_spectra       | str or dnpdata             | "all"       | "all"  unless dnpdata given                                          |
    +------------------+----------------------------+-------------+----------------------------------------------------------------------+
    | integrate_center | str, int, list, or ndarray | 0           | "max", center of integration window, or index used to find amplitude |
    +------------------+----------------------------+-------------+----------------------------------------------------------------------+
    | integrate_width  | str, int, list, or ndarray | "full"      | "full" or width of integration window                                |
    +------------------+----------------------------+-------------+----------------------------------------------------------------------+
    | method           | str                        | "integrate" | either "integrate" or "ampltiude"                                    |
    +------------------+----------------------------+-------------+----------------------------------------------------------------------+
    | dim              | str                        | "t2"        | dimension to integrate down or search down for max                   |
    +------------------+----------------------------+-------------+----------------------------------------------------------------------+

    Returns:
        all_data (dnpdata, dict): data object with "enhancement" added to the workspace

    """

    orig_data, isDict = return_data(all_data)
    if dim == "t2":
        ind_dim = "t1"
    elif dim == "t1":
        ind_dim = "t2"

    if (
        isinstance(off_spectrum, dnpdata)
        or isinstance(off_spectrum, dnpdata_collection)
    ) and (
        isinstance(on_spectra, dnpdata) or isinstance(on_spectra, dnpdata_collection)
    ):

        data_off, is_ws_off = return_data(off_spectrum)
        data_on, is_ws_on = return_data(on_spectra)

        if integrate_width == "full":
            int_width_off = "full"
            int_width_on = "full"
        else:
            if (
                isinstance(integrate_width, list)
                or isinstance(integrate_width, np.ndarray)
            ) and len(integrate_width) == 2:
                int_width_off = integrate_width[0]
                int_width_on = integrate_width[1]
            elif isinstance(integrate_width, int):
                int_width_off = integrate_width
                int_width_on = integrate_width
            else:
                raise ValueError(
                    "integrate width must be an integer, or a list/array with len = 2"
                )

        if (
            isinstance(integrate_center, list)
            or isinstance(integrate_center, np.ndarray)
        ) and len(integrate_center) == 2:
            int_center_off = integrate_center[0]
            int_center_on = integrate_center[1]
        elif isinstance(integrate_center, int):
            int_center_off = integrate_center
            int_center_on = integrate_center
        else:
            raise ValueError(
                "integrate center must be an integer, or a list/array with len = 2"
            )

        if method == "integrate":
            off_data1 = dnpNMR.integrate(
                data_off,
                dim=dim,
                integrate_center=int_center_off,
                integrate_width=int_width_off,
            )
            off_data = off_data1.values

            on_data1 = dnpNMR.integrate(
                data_on,
                dim=dim,
                integrate_center=int_center_on,
                integrate_width=int_width_on,
            )
            on_data = on_data1.values

        elif method == "amplitude":
            on_data = []
            if integrate_center == "max":
                off_data = data_off.values[np.argmax(abs(data_off.values))]
                if len(data_on.shape) == 1:
                    on_data.append(data_on.values[np.argmax(abs(data_on.values))])
                else:
                    for indx in range(data_on.shape[-1]):
                        on_data.append(
                            data_on.values[np.argmax(abs(data_on.values[indx])), indx]
                        )
            else:
                off_data = data_off.values[int_center_off]
                if len(data_on.shape) == 1:
                    on_data.append(data_on.values[int_center_on])
                else:
                    for indx in range(data_on.shape[-1]):
                        on_data.append(data_on.values[int_center_on, indx])

        if data_on.ndim == 2:
            enh_coords_on = data_on.coords[ind_dim]
        else:
            enh_coords_on = np.array(range(data_on.shape[-1]))

    elif isinstance(off_spectrum, int) and on_spectra == "all":

        if orig_data.ndim == 1:
            raise ValueError("data is 1D, enhancement will be equal to 1 !!")

        if off_spectrum == 0:
            off_spectrum = 1

        if method == "integrate":
            if integrate_width == "full":
                int_width = "full"
            else:
                int_width = integrate_width

            dnpNMR.integrate(
                all_data,
                dim=dim,
                integrate_center=integrate_center,
                integrate_width=int_width,
            )
            data, _isDict = return_data(all_data)
            data_1 = data.values

        elif method == "amplitude":
            data_1 = []
            if integrate_center == "max":
                for indx in range(orig_data.shape[-1]):
                    data_1.append(
                        orig_data.values[np.argmax(abs(orig_data.values[indx])), indx]
                    )
            else:
                for indx in range(orig_data.shape[-1]):
                    data_1.append(orig_data.values[integrate_center, indx])

        off_data = data_1[off_spectrum - 1]
        if off_spectrum == 1:
            on_data = data_1[1:]
            enh_coords_on = orig_data.coords[ind_dim][1:]
        elif off_spectrum > 1:
            on_data_1 = data_1[: off_spectrum - 1]
            on_coords_1 = orig_data.coords[ind_dim][: off_spectrum - 1]
            on_data_2 = data_1[off_spectrum:]
            on_coords_2 = orig_data.coords[ind_dim][off_spectrum:]
            on_data = np.concatenate((on_data_1, on_data_2))
            enh_coords_on = np.concatenate((on_coords_1, on_coords_2))

    enh = np.real(np.array(on_data) / np.array(off_data))

    enhancementData = dnpdata(enh, [enh_coords_on], [ind_dim])

    if isDict:
        all_data["enhancement"] = enhancementData
        return all_data
    else:
        return enhancementData
