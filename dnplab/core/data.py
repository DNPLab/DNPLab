"""
DNPData object for storing N-dimensional data with coordinates
"""

import numpy as _np

from .base import ABCData
from ..version import __version__

version = __version__

core_attrs_list = ["nmr_frequency"]

_np.set_printoptions(threshold=15)


class DNPData(ABCData):
    """
    DNPData Class for handling dnp data

    The DNPData class is inspired by pyspecdata nddata object which handles n-dimensional data, axes, and other relevant information together.

    This class is designed to handle data and axes together so that performing NMR processing can be performed easily.

    Attributes:
        values (numpy.ndarray): Numpy Array containing data
        coords (list): List of numpy arrays containing axes of data
        dims (list): List of axes labels for data
        attrs (dict): Dictionary of parameters for data


    """

    def __init__(
        self, values=_np.r_[[]], dims=[], coords=[], attrs = {}, dnplab_attrs = {}, proc_attrs={}
    ):
        """
        DNPData Class __init__ method

        Args:
            data (numpy.ndarray):
            coords (list): list of axes
            dims (list): list of strings which are names of axes
            attrs (dict): dictionary of parameters
            exp_attrs (dict): dictionary of experiment parameters
            dnplab_attrs (dict): dictionary of parameters used in dnplab 
            pro_attrs (dict): dictionary of parameters used in data processing
        """

        if len(dims) > 0 and isinstance(dims[0], list):
            dims = dims[0]

        super().__init__(values, dims, coords, attrs, dnplab_attrs)
        self.version = version
        self.attrs = attrs
        self.exp_attrs = attrs
        self.dnplab_attrs = dnplab_attrs
        self.proc_attrs = proc_attrs
        self.max_print_attrs = 5
        self.print_values = False

    @property
    def _constructor(self):
        return DNPData

    def __repr__(self):
        """
        Representation of DNPData object
        """
        return "nddata(values = {}, coords = {}, dims = {}, attrs = {})".format(
            repr(self.values), repr(self.coords), repr(self.dims), repr(self.attrs)
        )

    def __str__(self):
        """
        String representation of DNPData object
        """

        string = "values:\n\t"
        string += " x ".join(map(str, self.shape))

        string += " {} ({})\n".format(type(self.values).__name__, self.values.dtype)

        if self.print_values is True:
            string += str(self.values) + "\n"

        string += "dims:\n\t"

        string += "{}\n".format(self.dims)

        string += "coords:\n\t"
        string += "\n\t".join(map(repr, self.coords))

        string += "\n"

        string += "attrs:\n"

        for ix, key in enumerate(self.attrs.keys()):
            if ix == self.max_print_attrs:
                string += "\t+%i attrs" % (len(self.attrs) - self.max_print_attrs)
                break
            string += "\t{!r}: {!r}\n".format(key, self.attrs[key])

        return string

    def proc_info(self):
        """
        Print processing steps and parameters currently in proc_attrs list
        """

        print("-----------------")
        print("Processing Attributes")
        print("-----------------")
        if not self.proc_attrs:
            print("none.")
        else:
            longest_key = max(self.proc_attrs, key = len)
            maximum_length_of_proc_attrs = len(longest_key)
            for x in self.proc_attrs:
                spaces = " " * (1 + maximum_length_of_proc_attrs - len(x))
                print(
                    '| '
                    + x
                    + spaces
                    + '| '
                    + str(self.proc_attrs[x])
                    .replace("{", "")
                    .replace("}", "")
                    .replace("'", "")
                    )

    def exp_info(self):
        """
        Print experiment attributes currently in attrs dictionary
        """
        print("-----------------")
        print("Experiment Attributes")
        print("-----------------")
        if not self.attrs:
            print("none.")
        else:
            longest_key = max(self.attrs, key = len)
            maximum_length_of_attrs = len(longest_key)
            for x in self.attrs:
                spaces = " " * (1 + maximum_length_of_attrs - len(x))
                print(
                    '| '
                    + x
                    + spaces
                    + '| '
                    + str(self.attrs[x])
                    )
    
    def show_attrs(self, show_exp_info = True, show_proc_info = True):
        """
        Print experiment attributes and processing steps
        """

        if show_exp_info:
            self.exp_info()
        
        if show_proc_info:
            self.proc_info()

    def add_proc_attrs(self, proc_attr_name, proc_dict):
        """
        Stamp processing step to DNPData object

        Args:
            proc_attr_name (str): Name of processing step (e.g. "fourier_transform")
            proc_dict (dict): Dictionary of processing parameters for this processing step.
        """
        if not isinstance(proc_attr_name, str):
            raise ValueError("Processing step name must be string")
        if not isinstance(proc_dict, dict):
            raise ValueError("Processing dictionary must be dictionary")

        # self.proc_attrs.append((proc_attr_name, proc_dict))
        self.proc_attrs[proc_attr_name] = proc_dict

    def phase(self):
        """
        Return phase of DNPData object

        Returns:
            phase (float,int): phase of data calculated from sum of imaginary
                divided by sum of real components
        """
        return _np.arctan(
            _np.sum(_np.imag(self.values)) / _np.sum(_np.real(self.values))
        )

    def squeeze(self):
        """
        Remove all length 1 dimensions from data

        .. warning::
            Axes information is lost

        Example:
            data.squeeze()
        """
        remove_axes = []
        for axes_ix, axes_value in enumerate(self.coords):
            if len(axes_value) == 1:
                remove_axes.append(axes_ix)

        reverse_remove_axes = remove_axes[::-1]
        for index_ix, index_value in enumerate(reverse_remove_axes):
            self.coords.pop(index_value)
            self.dims.pop(index_value)
            self.values = _np.squeeze(self.values)

    def select(self, selection):
        """
        Select subset of 2D data object

        Args:
            selection (int, range, list, tuple): list or tuple of slices to keep

        Returns:
            DNPData object: subset of DNPData object

        Example:
            data.select((1, range(5,10), 15)) # keeps slices: 1, 5, 6, 7, 8, 9, and 15

        """
        if len(self.dims) == 1:
            raise TypeError("select method is not applicable to 1D data")
        if isinstance(selection, int):
            self.values = self.values[:, selection]
            self.coords[self.dims[1]] = self.coords[self.dims[1]][selection]
            self.dims = self.dims[0]
        elif isinstance(selection, range):
            self.values = self.values[:, selection.start : selection.stop]
            self.coords[self.dims[1]] = self.coords[self.dims[1]][
                selection.start : selection.stop
            ]
        elif isinstance(selection, (list, tuple)) and all(
            [isinstance(x, (int, range)) for x in selection]
        ):
            new_values = _np.empty(shape=(self.shape[0],))
            new_coords = []
            for x in selection:
                if isinstance(x, int):
                    new_values = _np.vstack((new_values, self.values[:, x]))
                    new_coords = new_coords + [self.coords[self.dims[1]][x]]
                elif isinstance(x, range):
                    new_values = _np.vstack(
                        (new_values, self.values[:, x.start : x.stop].T)
                    )
                    new_coords = new_coords + [
                        self.coords[self.dims[1]][y] for y in range(x.start, x.stop)
                    ]

            self.values = new_values.T[:, 1:]
            self.coords[self.dims[1]] = _np.array(new_coords)
        else:
            raise TypeError(
                "Select using integer, range, or list/tuple of integers or ranges"
            )

        return self._constructor(
            values=self.values,
            coords=self.coords._coords,
            dims=self.dims,
            attrs=self.attrs,
            proc_attrs=self.proc_attrs,
        )
