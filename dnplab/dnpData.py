"""dnpdata object for storing N-dimensional data with coordinates
"""
import numpy as np
from collections.abc import MutableMapping
from copy import deepcopy

from .core import nddata

version = "1.0"

core_attrs_list = ["nmr_frequency"]


class dnpdata(nddata.nddata_core):
    """dnpdata Class for handling dnp data

    The dnpdata class is inspired by pyspecdata nddata object which handles n-dimensional data, axes, and other relevant information together.

    This class is designed to handle data and axes together so that performing NMR processing can be performed easily.

    Attributes:
    values (numpy.ndarray): Numpy Array containing data
    coords (list): List of numpy arrays containing axes of data
    dims (list): List of axes labels for data
    attrs (dict): Dictionary of parameters for data

    """

    def __init__(self, values=np.r_[[]], coords=[], dims=[], attrs={}, procList=[]):
        """dnpdata Class __init__ method

        Args:
            data (numpy.ndarray):
            coords (list): list of axes
            dims (list): list of strings which are names of axes
            attrs (dict): dictionary of parameters


        """
        super().__init__(values, dims, coords, attrs)
        self.version = version
        self.proc_attrs = []

    #
    def __repr__(self):
        """Representation of dnpdata object"""
        return "nddata(values = {}, coords = {}, dims = {}, attrs = {})".format(
            repr(self.values), repr(self.coords), repr(self.dims), repr(self.attrs)
        )

    def __str__(self):
        """String representation of dnpdata object"""
        if len(self.attrs) < 20:
            return "values:\n{}\ndims:\n{}\ncoords:\n{}\nattrs:\n{}\nproc_attrs:\n{}".format(
                self.values, self.dims, self.coords, self.attrs, self.proc_attrs
            )
        else:
            core_attrs = {
                k: self.attrs[k] for k in core_attrs_list if k in core_attrs_list
            }
            num_additional_attrs = len(self.attrs) - len(core_attrs)
            return (
                "values:\n{}\ndims:\n{}\ncoords:\n{}\nattrs:\n{}\n + {} attrs".format(
                    self.values,
                    self.dims,
                    self.coords,
                    core_attrs,
                    num_additional_attrs,
                )
            )

    def add_proc_attrs(self, proc_attr_name, proc_dict):
        """Stamp processing step to dnpdata object

        Args:
            proc_attr_name (str): Name of processing step (e.g. "fourier_transform"
            proc_dict (dict): Dictionary of processing parameters for this processing step.
        """
        if not isinstance(proc_attr_name, str):
            raise ValueError("Processing step name must be string")
        if not isinstance(proc_dict, dict):
            raise ValueError("Processing dictionary must be dictionary")

        self.proc_attrs.append((proc_attr_name, proc_dict))

    def autophase(self):
        """Multiply dnpdata object by phase"""
        p = self.phase()
        self.values *= np.exp(-1j * p)
        if np.sum(np.real(self.values)) < 0:
            self.values *= -1.0

    def phase(self):
        """Return phase of dnpdata object

        Returns:
            phase (float,int): phase of data calculated from sum of imaginary divided by sum of real components
        """
        return np.arctan(np.sum(np.imag(self.values)) / np.sum(np.real(self.values)))

    def squeeze(self):
        """Remove all length 1 dimensions from data

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
            self.values = np.squeeze(self.values)


class dnpdata_collection(MutableMapping):
    """Dictionary-like workspace object for storing dnpdata objects"""

    def __init__(self, *args, **kwargs):
        """dnpdata_collection __init__ method

        Args:

        Example::

        """
        self._processing_buffer = "proc"

        self.__data_dict = {}

        if len(args) == 0:
            return
        elif len(args) == 1:
            if isinstance(args[0], dnpdata):
                self.__data_dict["raw"] == dnpdata
            elif isinstance(args[0], dict):
                data_dict = args[0]
                for key in data_dict:
                    if isinstance(data_dict[key], (dnpdata, dict)):
                        self.__data_dict[key] = data_dict[key]
                    else:
                        raise TypeError("Each type in dict must be dnpdata or dict")
            else:
                raise TypeError("Argument must be type dnpdata")
        elif len(args) == 2:
            if isinstance(args[0], str) and isinstance(args[1], (dnpdata, dict)):
                self.__data_dict[args[0]] = args[1]
            else:
                raise TypeError(
                    "If two arguments, first argument must be str and 2nd argument must be dnpdata or dict"
                )
        else:
            raise TypeError("Arguments not understood")

    def __getitem__(self, key):
        return self.__data_dict[key]

    def __setitem__(self, key, value):
        if (not isinstance(key, str)) or (not isinstance(value, (dict, dnpdata))):
            raise TypeError("Key must be string and value must be dnpdata or dict")
        self.__data_dict[key] = value

    def __delitem__(self, key):
        del self.__data_dict[key]

    def __iter__(self):
        return iter(self.__data_dict)

    def __len__(self):
        return len(self.__data_dict)

    @property
    def processing_buffer(self):
        return self._processing_buffer

    @processing_buffer.setter
    def processing_buffer(self, new_processing_buffer):
        """"""
        if isinstance(new_processing_buffer, str):
            self._processing_buffer = new_processing_buffer
        else:
            raise TypeError(
                "Processing buffer must be type str, not %s"
                % str(type(new_processing_buffer))
            )

    def copy(self, key, new_key=None):
        """Copy data from key to new_key. If new_key is not given, by default key will be copied to processing buffer

        Args:
            key (str): Key to be copied
            new_key (str, None): New key for copied data

        """

        if new_key is None:
            new_key = self.processing_buffer

        self[new_key] = self[key].copy()

    def move(self, key, new_key):
        """Move data from key to new_key

        Args:
            key (str): Name of data to move
            new_key (str): Name of new key to move data
        """

        self[new_key] = self.pop(key)

    def pop(self, key):
        """Pop key. Removes data corresponding to key."""
        return self.__data_dict.pop(key)

    def dict(self):
        """Return dictionary for storing data in dnpdata_collection"""
        return self.__data_dict

    def clear(self):
        """Removes all items"""
        self.__data_dict.clear()

    get = __getitem__

    def items(self):
        """Return items"""
        return self.__data_dict.items()

    def keys(self):
        """Return keys."""
        return self.__data_dict.keys()

    def popitem(self):
        """Pops item from end of dnpdata_collection

        Returns:
            tuple: key, item pair that was removed
        """
        return self.__data_dict.popitem()

    def values(self):
        """Return Values"""
        return self.__data_dict.values()

    def add(self, key, data):
        """Adds new data

        Args:
            key (str): key corresponding to new data
            data (dnpdata): data object corresponding to key
        """
        if (not isinstance(key, str)) or (not isinstance(data, (dnpdata, dict))):
            raise TypeError("add takes two arguments, a string and dnplab.dnpdata type")
        self.__data_dict[key] = data

    def __repr__(self):
        return "dnpdata_collection({})".format(self.__data_dict)

    def __str__(self):
        return "{}\n".format([(key, self[key].__str__()) for key in self.keys()])


def create_workspace(*args):
    """Create a workspace (dnpdata_collection)

    Args:
        args: Arguments to send to __init__ method in dnpdata_collection

    Returns:
        dnpdata_collection: workspace object
    """
    return dnpdata_collection(*args)


if __name__ == "__main__":
    pass
