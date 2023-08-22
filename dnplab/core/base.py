from __future__ import division
import operator
import numpy as _np
import warnings
from copy import deepcopy
from collections import OrderedDict
from .coord import Coords
import logging
import warnings

from ..version import __version__

version = __version__

logger = logging.getLogger(__name__)


def _replaceClassWithAttribute(replace_class, args, kwargs, target_attr="_values"):
    r_args = []
    for a in args:
        if type(a) == type(replace_class):
            r_args.append(getattr(replace_class, target_attr))
        else:
            r_args.append(a)
    r_kwargs = {}
    # iterating over .values() can lead to errors?
    for key in kwargs.keys():
        val = kwargs[key]
        if type(val) == type(replace_class):
            r_kwargs[key] = getattr(replace_class, target_attr)
        else:
            r_kwargs[key] = val
    return tuple(r_args), r_kwargs


# utility funtion to return integer index of string or int if input is integer
_str_to_int_index = (
    lambda possible_dim, dnpdat: int(dnpdat.index(possible_dim))
    if isinstance(possible_dim, str)
    else possible_dim
)

_SPECIAL_NP_HANDLED = {}  # numpy functions that need special handling, will

_numerical_types = (
    _np.ndarray,
    int,
    float,
    complex,
    _np.complex64,
    _np.float16,
    _np.float32,
    _np.double,
    _np.longdouble,
    _np.csingle,
    _np.cdouble,
    _np.clongdouble,
    _np.longlong,
)

# _nddata_core_version = "1.0"


class ABCData(object):
    """N-Dimensional Data Object

    Attributes:
        values (numpy.ndarray): Data values in
        dims (list): List of strings giving dimension labels
        coords (Coords): Collection of numpy.ndarrays defining the axes
        attrs (dict): dictionary of parameters
        error (numpy.ndarray): If not None, error for values which are propagated during mathematical operations
        proc_attrs (list): List of processing steps
    """

    __array_priority__ = 1000  # radd, rsub, ... should return nddata object

    def __init__(
        self,
        values=_np.r_[[]],
        dims=[],
        coords=[],
        attrs={},
        dnp_attrs={},
        error=None,
        **kwargs
    ):
        self.version = version

        # if values is list, convert to numpy array
        if isinstance(values, list):
            values = _np.array(values)

        # verify values are numpy array
        if isinstance(values, _np.ndarray):
            self._values = values
        else:
            raise TypeError(
                'values must be type "numpy.ndarray" not %s' % str(type(values))
            )

        self._coords = Coords(dims, coords)

        if isinstance(attrs, dict):
            self._attrs = attrs
        else:
            raise TypeError('attrs must be type "dict" not %s' % str(type(attrs)))

        if isinstance(dnp_attrs, dict):
            self._attrs = attrs
        else:
            raise TypeError(
                'dnp_attrs must be type "dict" not %s' % str(type(dnp_attrs))
            )

        if isinstance(error, _np.ndarray) or (error == None):
            self._error = error
        else:
            raise TypeError(
                'error must be type "numpy.ndarray" or "None" not %s' % str(type(error))
            )

        if "proc_attrs" in kwargs:
            proc_attrs = kwargs["proc_attrs"]
            if isinstance(proc_attrs, list):
                self._proc_attrs = proc_attrs
            else:
                raise TypeError

        if not self._self_consistent():
            warnings.warn("Data Object Dimensions are not consistent.")
            if not self._check_dims(self.dims):
                warnings.warn(
                    "Dims Check Failed. Dims must be list of strings. Length of dims ({0}) should match number of dimensions in values (values shape: {1})".format(
                        len(self.dims), self._values.shape
                    )
                )
            if not self._check_coords(self._coords):
                warnings.warn(
                    "Coords Check Failed. Each coord should be a numpy array matching the length of each dimension in values."
                )
        self._is_folded = True

    @property
    def __version__(self):
        return self._nddata_core_version

    def _check_dims(self, dims):
        """Check that list is a list of strings with no duplicates

        Args:
            dims list to verify as list of strings

        Returns:
            bool: True if arguments is list of strings with no duplicates, False otherwise
        """

        # test that all members are strings
        all_strings = all([isinstance(dims[x], str) for x in range(len(dims))])

        # test if any duplicates exist
        any_duplicates = len(dims) == len(set(dims))

        # Test if number of dims matches number of dimensions
        if self.size != 0:  # Case where array is not empty
            shape_check = len(self.values.shape) == len(self.dims)
        else:  # Case where array is empty
            shape_check = len(self.dims) == 0

        return all_strings and any_duplicates and shape_check

    def _check_coords(self, coords):
        """Check that coords is list of 1d numpy arrays

        Args:
            coords: coords object to check if valid coords

        Returns:
            bool: True if valid coords. False, otherwise.
        """

        # Check that coords is list of 1d numpy arrays
        for coord in coords:
            if isinstance(coord, _np.ndarray):
                if len(coord.shape) == 1:
                    pass
                else:
                    return False
            else:
                return False

        shape = self.values.shape

        # Check Shapes are consistent
        if len(shape) != len(coords):
            return False

        # Check that each coord length matches shape of each dimension in values
        for coord_ix, coord in enumerate(coords):
            if len(coord) != shape[coord_ix]:
                return False

        return True

    def _check_error(self, error):
        """Check that error attribute is numpy array and it's size matches values

        Args:
            error: Object to check if valid error

        Returns:
            bool: True if valid error. False otherwise.
        """

        # Error is type None if it doesn't exist
        if error is None:
            return True

        check_type = isinstance(error, _np.ndarray)

        if check_type:
            check_size = error.size == self._values.size
        else:
            check_size = False

        return check_type and check_size

    def _self_consistent(self):
        """Check if dimensions for values, dims, and coords are consistent

        Returns:
            bool: True if consistent, False otherwise
        """

        if self._values.size == 0:
            coords_check = len(self._coords) == 0
        else:
            coords_check = list(self._values.shape) == list(self.coords.shape)

        dims_check = self._check_dims(self.dims)

        return coords_check and dims_check

    def _attrs_valid(self):
        """Verify attrs attribute is valid. All values in attrs must be list, numpy.ndarray, int, float, complex, or str.

        Returns:
            bool: True if attrs is valid. False, otherwise.
        """

        for key in self._attrs:
            if not isinstance(
                self._attrs, (list, _np.ndarray, int, float, complex, str)
            ):
                return False

        return True

    def __getitem__(self, args):
        """Method for indexing nddata_core

        Args:
            args (tuple): Tuple containing alternative dims and indexing values for each dimension to be indexed. (e.g. data['x', 1:10, 'y', :, 'z', (3.5, 7.5)])

        Example::

            data['x', 1] # return data indexing down "x" dim with index 1

            data['x', 4.5] # return data indexing down "x" dim with "x" coords nearest to 4.5

            data['x', 2:5] # return data indexing down "x" dim with index from 2 to 5

            data['x', (100., 150.)] # return data indexing down "x" dim where coords are range from 100 to 150
        """
        a = self.copy()
        if len(args) % 2 == 1:
            raise ValueError("Cannot index with odd number of arguments")

        index_dims = args[0::2]
        for dim in index_dims:
            if dim not in a.dims:
                raise ValueError("dim not in dims")

        index_slice = args[1::2]

        # check slices
        for slice_ in index_slice:
            # type must be slice or tuple
            if not isinstance(slice_, (slice, tuple, float, int)):
                raise ValueError("Invalid slice type")

            # if tuple, length must be two: (start, stop)
            if isinstance(slice_, tuple) and not len(slice_) in (1, 2):
                raise ValueError("tuple index must have one or two values")

        # convert tuple to slice
        updated_index_slice = []
        for dim, slice_ in zip(index_dims, index_slice):
            if isinstance(slice_, tuple):
                index = a.index(dim)
                if len(slice_) == 1:
                    start = _np.argmin(_np.abs(slice_[0] - a.get_coord(dim)))
                    updated_index_slice.append(slice(start, start + 1))
                else:
                    start = _np.argmin(_np.abs(slice_[0] - a.get_coord(dim)))
                    stop = _np.argmin(_np.abs(slice_[1] - a.get_coord(dim)))
                    if start == stop:
                        stop = start + 1
                    if stop < start:
                        start, stop = stop, start
                    updated_index_slice.append(slice(start, stop))
            elif isinstance(slice_, int):
                start = slice_
                if slice_ != -1:
                    updated_index_slice.append(slice(start, start + 1))
                else:
                    updated_index_slice.append(slice(slice_, None))
            elif isinstance(slice_, float):
                start = _np.argmin(_np.abs(slice_ - a.get_coord(dim)))
                updated_index_slice.append(slice(start, start + 1))
            else:
                updated_index_slice.append(slice_)

        # special behavior for unfolded (internal use):
        if a._is_folded == False:
            if len(index_dims) > 2:
                raise IndexError(
                    "Unfolded state, can request at maximum 2 dimensions, you requested {0}".format(
                        index_dims
                    )
                )
            for k in index_dims:
                if (not k in ["fold_index",self.dims[0]]):
                    raise IndexError(
                        "Unfolded state, there are only two dimensions, you specified {0} but you can only request {1}/{2}".format(
                            index_dims,self.dims[0],'fold_index')
                        )
            # now some random slices have been selected, it is only useful to consider this as a new object which is considered as folded
            # pop coords except self.dims[0] and "fold_index"
            for d in self.dims[1:]:
                if d=="fold_index":
                    continue
                a.coords.pop(d)
            index_slice_dict = dict(zip(index_dims, updated_index_slice))
            new_slices = [slice(None) if (dim not in index_dims) else index_slice_dict[dim] for dim in a.dims]
            a.values = a.values[tuple(new_slices)]

            a.coords.pop('fold_index')
            a.coords['fi']=_np.arange(a.shape[1])
            a._is_folded=True
            return a

        index_slice_dict = dict(zip(index_dims, updated_index_slice))
        new_slices = [
            slice(None) if dim not in index_dims else index_slice_dict[dim]
            for dim in a.dims
        ]

        for ix, dim in enumerate(a.dims):
            a.coords[dim] = a.coords[dim][new_slices[ix]]

        a.values = a.values[tuple(new_slices)]

        return a

    def __setitem__(self, args, new_values):
        """Method for setting values by index nddata_core

        Args:
            args (tuple): Tuple containing alternative dims and indexing values for each dimension to be indexed. (e.g. data['x', 1:10, 'y', :, 'z', (3.5, 7.5)])
            new_values (numpy.ndarray): New values at given index

        Example::

            data['x', 1] = 1 # return data indexing down "x" dim with index 1

        """
        a = self.copy()

        if isinstance(new_values, ABCData):
            new_values = new_values.values

        if len(args) % 2 == 1:
            raise ValueError("Cannot index with odd number of arguments")

        index_dims = args[0::2]
        for dim in index_dims:
            if dim not in a.dims:
                raise ValueError("dim not in dims")

        index_slice = args[1::2]

        # check slices
        for slice_ in index_slice:
            # type must be slice or tuple
            if not isinstance(slice_, (slice, tuple, float, int)):
                raise ValueError("Invalid slice type")

            # if tuple, length must be two: (start, stop)
            if isinstance(slice_, tuple) and not len(slice_) in (1, 2):
                raise ValueError("tuple index must have one or two values")

        # convert tuple to slice
        updated_index_slice = []
        for dim, slice_ in zip(index_dims, index_slice):
            if isinstance(slice_, tuple):
                index = a.index(dim)
                if len(slice_) == 1:
                    start = _np.argmin(_np.abs(slice_[0] - a.get_coord(dim)))
                    updated_index_slice.append(slice(start, start + 1))
                else:
                    start = _np.argmin(_np.abs(slice_[0] - a.get_coord(dim)))
                    stop = _np.argmin(_np.abs(slice_[1] - a.get_coord(dim)))
                    if start == stop:
                        stop = start + 1
                    if stop < start:
                        start, stop = stop, start
                    updated_index_slice.append(slice(start, stop))
            elif isinstance(slice_, int):
                start = slice_
                if slice_ != -1:
                    updated_index_slice.append(slice(start, start + 1))
                else:
                    updated_index_slice.append(slice(slice_, None))
            elif isinstance(slice_, float):
                start = _np.argmin(_np.abs(slice_ - a.get_coord(dim)))
                updated_index_slice.append(slice(start, start + 1))
            else:
                updated_index_slice.append(slice_)

        index_slice_dict = dict(zip(index_dims, updated_index_slice))
        new_slices = [
            slice(None) if dim not in index_dims else index_slice_dict[dim]
            for dim in a.dims
        ]

        self.values[tuple(new_slices)] = new_values

    def copy(self):
        """Return deepcopy of dnpdata object

        Returns:
            deep copy of data object
        """
        return deepcopy(self)

    def cumulative_sum(self, dim):
        """Calculate Cumulative sum of dnpdata object

        Returns:
            cumulative sum of data object
        """
        a = self.copy()
        index = a.index(dim)
        a.values = a.values.cumsum(index)
        # NOTE: Add Error Propagation
        return a

    def merge_attrs(self, b):
        """Merge the given dictionaries

        Args:
            b (nddata_core): attributes to merge into object
        """

        for key in b.attrs:
            if key not in self.attrs:
                self.attrs[key] = b.attrs[key]
            else:
                if not self.attrs[key] is b.attrs[key]:
                    pass
                    # warnings.warn(
                    #     "attrs in two dictionarys contain different values, leaving original value:\n{}:{}".format(
                    #         key, self.attrs[key]
                    #     )
                    # )

    def __len__(self):
        """Returns len(values) the length of the first dimension in values"""
        return len(self._values)

    @property
    def size(self):
        """Returns values.size. Total number of elements in numpy array."""
        return self._values.size

    def sort_dims(self):
        """Sort the dimensions"""
        sorted_order = sorted(range(len(self.dims)), key=lambda x: self.dims[x])

        self.coords.reorder_index(sorted_order)
        self._values = _np.moveaxis(
            self._values, range(len(sorted_order)), sorted_order
        )

    def index(self, dim):
        """Find index of given dimension name

        Args:
            dim (str): Name of dimension to index

        Returns:
            int: Index value of dim
        """
        if dim in self.dims:
            return self.coords.index(dim)
        else:
            raise ValueError(
                'The dimension "%s" is not in dims. Available dimensions are: %s'
                % (dim, self.dims)
            )

    def __truediv__(self, b):
        if isinstance(b, ABCData):
            a, b = self.align(b)

            a.values = a.values / b.values

            # error propagation
            if a.error is not None and b.error is not None:
                error = abs(a.values) * _np.sqrt(
                    (self.error / a.values) ** 2.0 + (b.error / a.values) ** 2.0
                )
            elif b.error is not None:
                a.error = b.error

            return a

        elif isinstance(b, _numerical_types):
            a = self.copy()
            a.values = a.values.__truediv__(b)
            return a
        else:
            raise TypeError("Cannot add type: {}".format(type(b)))

    def __rtruediv__(self, b):
        if isinstance(b, _numerical_types):
            a = self.copy()
            a.values = self.values.__rtruediv__(b)
            return a
        else:
            raise TypeError("Cannot add type: {}".format(type(b)))

    def __mul__(self, b):
        if isinstance(b, ABCData):
            a, b = self.align(b)

            a.values = a.values * b.values

            # error propagation
            if a.error is not None and b.error is not None:
                error = abs(a.values) * _np.sqrt(
                    (self.error / a.values) ** 2.0 + (b.error / a.values) ** 2.0
                )
            elif b.error is not None:
                a.error = b.error

            return a

        elif isinstance(b, _numerical_types):
            a = self.copy()
            a.values = a.values.__mul__(b)
            return a
        else:
            raise TypeError("Cannot add type: {}".format(type(b)))

    __rmul__ = __mul__

    def __add__(self, b):
        if isinstance(b, ABCData):
            a, b = self.align(b)

            a.values = a.values + b.values

            # error propagation
            if a.error is not None and b.error is not None:
                a.error = _np.sqrt(a.error**2.0 + b.error**2.0)
            elif b.error is not None:
                a.error = b.error

            return a

        elif isinstance(b, _numerical_types):
            a = self.copy()
            a.values = a.values.__add__(b)
            return a
        else:
            raise TypeError("Cannot add type: {}".format(type(b)))

    __radd__ = __add__

    def __sub__(self, b):
        if isinstance(b, ABCData):
            a, b = self.align(b)

            a.values = a.values - b.values

            # error propagation
            if a.error is not None and b.error is not None:
                a.error = _np.sqrt(a.error**2.0 + b.error**2.0)
            elif b.error is not None:
                a.error = b.error

            return a

        elif isinstance(b, _numerical_types):
            a = self.copy()
            a.values = a.values.__sub__(b)
            return a
        else:
            raise TypeError("Cannot add type: {}".format(type(b)))

    def __rsub__(self, b):
        if isinstance(b, _numerical_types):
            a = self.copy()
            a.values = a.values.__rsub__(b)
            return a
        else:
            raise TypeError("Cannot add type: {}".format(type(b)))

    @property
    def values(self):
        return self._values

    @values.setter
    def values(self, b):
        if not isinstance(b, (int, complex, float, _np.ndarray)):
            raise TypeError('Values must be type "numpy.ndarray" not %s' % type(b))
        if isinstance(b, (int, complex, float)):
            b = _np.array(b)

        self._values = b

    @values.getter
    def values(self):
        return self._values

    @property
    def dims(self):
        return self.coords.dims

    @dims.setter
    def dims(self, b):
        self.coords.dims = b

    @property
    def coords(self):
        return self._coords

    @coords.setter
    def coords(self, b):
        self._coords = b

    @property
    def attrs(self):
        return self._attrs

    @attrs.setter
    def attrs(self, b):
        if isinstance(b, dict):
            self._attrs = b
        else:
            raise ValueError('attrs must be type "dict"')

    @property
    def error(self):
        return self._error

    @error.setter
    def error(self, b):
        if isinstance(b, (_np.ndarray, type(None))):
            self._error = b
        else:
            raise ValueError('error must be type "numpy.ndarray"')

    def rename(self, dim, new_name):
        """Rename dim

        Args:
            dim (str): Name of dimension to rename
            new_name (str): New name for dim
        """
        if dim in self.dims:
            if isinstance(new_name, str):
                self.coords.rename(dim, new_name)
            else:
                raise TypeError(
                    'New dimension name must be type "str" not %s' % type(new_name)
                )
        else:
            raise ValueError("Dimension name %s is not in dims" % dim)

    def reorder(self, dims):
        """Reorder dimensions

        Args:
            dims (list): List of strings in new order

        """

        if not self._check_dims(dims):
            raise TypeError("New dims must be list of str with no duplicates")
        for dim in dims:
            if dim not in self.dims:
                raise ValueError("no such dimension: %s" % dim)

        # Add original dims to end, remove duplicates
        dims = list(OrderedDict.fromkeys(dims + self.dims))

        #        print('original dims', self.dims)
        #        print('new order of dims', dims)

        new_order = [dims.index(dim) for dim in self.dims]
        permutation_order = [self.dims.index(dim) for dim in dims]
        #        print('New Order', new_order)
        #        print('Permutation Order', permutation_order)

        self.coords.reorder(dims)
        #        print(self.coords)

        # Transpose values
        #        self.values = np.transpose(self.values, new_order)
        self.values = _np.transpose(self.values, permutation_order)

    def __str__(self):
        return "values:\n{}\ndims:\n{}\ncoords:\n{}\nattrs:\n{}".format(
            self.values, self.dims, self.coords, self.attrs
        )

    def __repr__(self):
        return "nddata_core(values = {}, coords = {}, dims = {}, attrs = {})".format(
            repr(self.values), repr(self.dims), repr(self.coords), repr(self.attrs)
        )

    def squeeze(self, dim):
        """Remove length 1 axes"""
        a = self.copy()
        shape = a.shape

        remove_dims = [a.dims[x] for x in range(len(shape)) if shape[x] == 1]
        values = _np.squeeze(a.values)

        if a.error is not None:
            a.error = _np.squeeze(a.error)

        attrs = a.attrs

        for dim in remove_dims:
            out = a.coords.pop(dim)
            if dim not in attrs:
                attrs[dim] = _np.array(out)
            else:
                warnings.warn("Attribute lost {}:{}".format(dim, out))

        return a

    def chunk(self, dim, new_dims, new_sizes):
        """
        .. note::
            This is a placeholder for a function that's not yet implemented

        Parameters
        ----------
        dim: str
            Assume that the dimension `dim` is a direct product of the
            dimensions given in `new_dims`, and chunk it out into those new
            dimensions.
        new_dims: list of str
            The new dimensions to generate.  Note that one of the elements of
            the list can be `dim` if you like.

            It's assumed that the ordering of `dim` is a direct product given
            in C-ordering (*i.e.* the inner dimensions are listed last and the
            outer dimensions are listed first -- here "inner" means that
            changes to the index of the inner-most dimension correspond to
            adjacent positions in memory and/or adjacent indeces in the
            original dimension that you are chunking)
        new_sizes: list of int
            sizes of the new dimensions`
        Returns
        -------
        self: nddata_core
            The new nddata object.
            Note that uniformly ascending or descending coordinates are manipulated in a rational way,
            *e.g.* `[1,2,3,4,5,6]` when chunked to a size of `[2,3]` will yield
            coordinates for the two new dimensions:
            `[1,4]` and `[0,1,2]`.
            Coordinates that are not uniformly ascending or descending will
            yield and error and must be manually modified by the user.
        """
        return NotImplemented

    def get_coord(self, dim):
        """Return coord corresponding to given dimension name

        Args:
            dim (str): Name of dim to retrieve coordinates from

        Returns:
            numpy.ndarray: array of coordinates

        """
        return self.coords[dim]

    @property
    def shape(self):
        """tuple: Shape of values"""
        return self.values.shape

    @property
    def dtype(self):
        """type: Values type"""
        return self.values.dtype

    def sum(self, dim):
        """Perform sum down given dimension

        Args:
            dim (str): Dimension to perform sum down

        """
        a = self.copy()
        index = a.index(dim)
        a.values = a.values.sum(index)
        if a.error is not None:
            a.error = a.error.std(index)
        a.coords.pop(dim)
        DeprecationWarning("Method '.sum()' will be removed on September 1st, 2023 ")
        return a

    def align(self, b):
        """Align two data objects for numerical operations

        Args:
            b: Object to align with self

        Returns:
            tuple: self and b aligned data objects
        """
        a = self.copy()
        b = b.copy()

        all_dims = list(OrderedDict.fromkeys(a.dims + b.dims))
        new_b_order = [dim for dim in all_dims if dim in b.dims]
        new_order = [b.index(dim) for dim in new_b_order]

        # create new dims where necessary
        values = a.values[
            tuple(slice(None) if dim in a.dims else None for dim in all_dims)
        ]

        # re-order
        old_order = list(range(len(new_order)))
        # re-order b values so they match order of all_dims
        values_b = _np.moveaxis(b.values, new_order, old_order)
        # create new dims where necessary
        values_b = values_b[
            tuple(slice(None) if dim in b.dims else None for dim in all_dims)
        ]

        # Handle Error
        if a.error is not None:
            error = a.error[
                tuple(slice(None) if dim in a.dims else None for dim in all_dims)
            ]
        else:
            error = a.error
        if b.error is not None:
            error_b = _np.moveaxis(b.values, new_order, old_order)
            error_b = error_b[
                tuple(slice(None) if dim in b.dims else None for dim in all_dims)
            ]
        else:
            error_b = b.error

        # check coords
        for dim in all_dims:
            if (dim in a.dims) and (dim in b.dims):
                if not _np.allclose(a.get_coord(dim), b.get_coord(dim)):
                    raise ValueError("Coords do not match for dim: %s" % dim)

        # merge attrs
        a.merge_attrs(b)

        # coords
        coords = a.coords.copy()
        for dim in new_b_order:
            if dim not in a.dims:
                coords.append(dim, b.coords[dim])

        a.values = values
        a.coords = coords
        a.error = error

        b.values = values_b
        b.coords = coords
        b.error = error_b

        return a, b

    def __array__(self):
        return self.values

    def smoosh(self, old_dims, new_name):
        """
        .. note::
            Not yet implemented.

        `smoosh` does the opposite of `chunk` -- see :func`:~nddata_core.chunk`
        """
        return NotImplemented

    def sort(self, dim):
        """Sort the coords corresponding to the given dim in ascending order

        Args:
            dim (str): dimension to sort
        """

        sort_array = _np.argsort(self.coords[dim])

        self.coords[dim] = self.coords[dim][sort_array]

        new_order = tuple(
            [slice(None) if dim != this_dim else sort_array for this_dim in self.dims]
        )

        self.values = self.values[new_order]

    def split(self, dim, new_dim, coord):
        """Split the dimension dim into"""

        if isinstance(coord, int):
            coord = _np.arange(coord)

        # move dim to end of dims
        dims = self.dims
        dims.remove(dim)
        dims.append(dim)

        # reorder data with split dim at end
        self.reorder(dims)

        new_shape = list(self.coords.shape)
        new_shape[-1] = int(new_shape[-1] / coord.size)
        self.coords[dim] = self.coords[dim][0 : new_shape[-1]]

        new_shape += [coord.size]
        print(new_shape)

        self.values = self.values.reshape(new_shape)
        self.coords.append(new_dim, coord)

    def is_sorted(self, dim):
        """Determine if coords corresponding to give dim are sorted in ascending order
        Args:
            dim (str): Dimension to check if sorted

        Returns:
            bool: True if sorted, False otherwise.
        """
        return _np.all(self.coords[dim][:-1] <= self.coords[dim][1:])

    @property
    def real(self):
        """DNPData: DNPData with real part of values"""
        a = self.copy()
        a.values = _np.real(a.values)
        return a

    @property
    def imag(self):
        """DNPData: DNPData with imaginary part of values"""
        a = self.copy()
        a.values = _np.imag(a.values)
        return a

    @property
    def abs(self):
        """DNPData: DNPData with absolute part of values"""
        a = self.copy()
        a.values = _np.abs(a.values)
        return a

    def concatenate(self, b, dim):
        """Concatenate DNPData objects

        Args:
            b (DNPData): Data object to append to current data object
            dim (str): dimension to concatenate along
        """

        if not dim in b.dims:
            raise ValueError("dim does not exist")
        if not dim in self.dims:
            raise ValueError("dim does not exist")

        index = self.dims.index(dim)

        b.reorder(self.dims)

        self.values = _np.concatenate((self.values, b.values), axis=index)

        self.coords[dim] = _np.concatenate(
            (
                _np.array(self.coords[dim]).reshape(-1),
                _np.array(b.coords[dim]).reshape(-1),
            )
        )

    def new_dim(self, dim, coord):
        """Add new dimension with length 1

        Args:
            dim (str): Name of new dimension
            coord (int, float): New coord
        """
        self.coords.append(dim, _np.r_[coord])
        self.values = _np.expand_dims(self.values, -1)

    def maximum(self, dim):
        """Return max for given dim

        Args:
            dim (str): Dimension to take maximum along

        """
        a = self.copy()
        index = a.dims.index(dim)

        a.values = _np.max(a.values, axis=index)
        a.coords.pop(dim)

        return a

    def argmax(self, dim):
        """Return value of coord at values maximum for given dim

        Args:
            dim (str): Dimension to perform operation along
        """
        a = self.copy()
        index = a.dims.index(dim)

        a.values = a.coords[dim][_np.argmax(a.values, axis=index)]
        a.coords.pop(dim)

        return a

    def argmax_index(self, dim):
        """Return index of coord at values maximum for given dim

        Args:
            dim (str): Dimension to perform operation along
        """
        a = self.copy()
        index = a.dims.index(dim)

        a.values = _np.argmax(a.values, axis=index)
        a.coords.pop(dim)

        return a

    def minimum(self, dim):
        """Return min for given dim

        Args:
            dim (str): Dimension to perform operation along
        """
        a = self.copy()
        index = a.dims.index(dim)

        a.values = _np.min(a.values, axis=index)
        a.coords.pop(dim)

        return a

    def argmin(self, dim):
        """Return value of coord at values minimum for given dim

        Args:
            dim (str): Dimension to perform operation along
        """
        a = self.copy()
        index = a.dims.index(dim)

        a.values = a.coords[dim][_np.argmin(a.values, axis=index)]
        a.coords.pop(dim)

        return a

    def argmin_index(self, dim):
        """Return index of coord at values minimum for given dim

        Args:
            dim (str): Dimension to perform operation along
        """
        a = self.copy()
        index = a.dims.index(dim)

        a.values = _np.argmin(a.values, axis=index)
        a.coords.pop(dim)

        return a

    @property
    def ndim(self):
        """str: Number of dimensions"""
        return self.values.ndim

    def unfold(self, dim):
        """Unfold ND data to 2d data

        Args:
            dim (str): Dimension to make first (length N), all other dimensions unfolded so that values has shape (N x M)

        """
        if self._is_folded == False:
            raise ValueError(
                "Trying to unfold already unfolded object (self._is_folded == False)! Please fold the object before attempting a second unfold!"
            )

        folded_order = self.dims  # Original order of dims
        self.reorder([dim])  # Move dim to first dimension
        folded_shape = self.values.shape
        align_dim_length = folded_shape[0]  # length of dimension to align down
        self.values = self.values.reshape(align_dim_length, -1)  # Reshape to 2d

        # remark: this will fail when unfolding an already unfolded dataset!
        index_dim = self.values.shape[1]
        self.coords.append("fold_index", _np.arange(index_dim))

        self.attrs["folded_shape"] = folded_shape
        self.attrs["folded_order"] = folded_order
        self._is_folded = False

    def fold(self):
        """Fold 2d data to original ND shape"""
        # remark: this will fail when unfolding an already unfolded dataset!
        self.coords.pop("fold_index")
        if "folded_shape" in self.attrs:
            original_shape = self.attrs["folded_shape"]
            self.values = self.values.reshape(original_shape)
            self.attrs.pop("folded_shape")
        if "folded_order" in self.attrs:
            self.reorder(self.attrs["folded_order"])
            self.attrs.pop("folded_order")
        self._is_folded = True

    def __array_ufunc__(self, ufunc, method, *arrInput, **kwargs):
        """
        started implementation of numpy comapilibity according to
        https://numpy.org/doc/stable/user/basics.dispatch.html

        - still need to implement universal function handling:
        https://numpy.org/doc/stable/user/basics.ufuncs.html#ufuncs-basics

        - need to implement all methods
        """
        if method == "__call__":
            args, kwargs = _replaceClassWithAttribute(self, arrInput, kwargs)
            values = ufunc(*args, **kwargs)
            if kwargs.pop("_dnplab_inplace", False):
                a = self
            else:
                a = self.copy()
            a.values = values
            proc_attr_name = "numpy." + ufunc.__name__
            proc_parameters = {"args": args, "kwargs": kwargs}
            try:
                a.add_proc_attrs(proc_attr_name, proc_parameters)
            except AttributeError:
                pass  # would log error here
            return a
        return NotImplemented

    def __array_function__(self, func, types, args, kwargs):
        """
        Basic numpy compability according to
        https://numpy.org/doc/stable/user/basics.dispatch.html

        - special handling is done for function that return no array, be aware that this is currently not optimal and possibly slow
        - numpy functions currently just apply to the values of self, special handling is done for functions in SPECIAL_NP_HANDLED

        - default implementation, replaces all  ABCData objects with ABCData.values and calls func with replaced *args,**kwargs
        - todo: better handling of axis keyword and dimension consume!
        """
        if func in _SPECIAL_NP_HANDLED:
            return_values = _SPECIAL_NP_HANDLED[func](*args, **kwargs)
        else:
            # default implementation, replaces all  ABCData objects in args and kwargs with ABCData.values and calls func
            args, kwargs = _replaceClassWithAttribute(self, args, kwargs)

            data_dims = []
            proc_dims = []
            # the following construction is needed as axis default value is not always None in numpy
            # handle axis keyword, if existent
            if "axis" in kwargs.keys():
                # pop axis keyword and if it is None stop here
                ax_value = kwargs.pop("axis")
                proc_dims.append(ax_value)
                if ax_value is None:
                    kwargs["axis"] = None
                elif type(ax_value) == int:
                    kwargs["axis"] = ax_value

                else:
                    # axis could now be a string or a tuple
                    if isinstance(ax_value, str):
                        # in case of string: find corresponding axis index and add to data dims
                        indx = tuple([int(self.index(ax_value))])
                        data_dims += [ax_value]
                        kwargs["axis"] = indx
                    else:
                        indx = []
                        for value in ax_value:
                            _val = _str_to_int_index(value)
                            if isinstance(value, str):
                                data_dims += [value]
                            indx += [_val]
                        indx = tuple(indx)
                        kwargs["axis"] = indx
            else:
                # use default value for axis argument
                pass
            # forbid out keyword as this is not implemented
            if "out" in kwargs.keys():
                raise NotImplemented

            # apply function to values
            return_values = func(*args, **kwargs)

        if type(return_values) == _np.ndarray:
            self_shape = self._values.shape
            if kwargs.pop("_dnplab_inplace", False):
                a = self
            else:
                a = self.copy()
            # delete dimensions that are removed
            # beware: case not handled: when only partial dimensions are deleted
            # current temporary fix: when shape of return_values.shape == self._values.shape do not delete dimensions!
            if return_values.shape != self_shape:
                for dim in data_dims:
                    a.coords.pop(dim)
            a.values = return_values

            # add processing attributes, except when they should be suppressed
            if kwargs.pop("_dnplab_suppress_attr", False):
                return a
            proc_attr_name = "numpy." + func.__name__
            proc_parameters = {"axis": proc_dims}
            try:
                a.add_proc_attrs(proc_attr_name, proc_parameters)
            except AttributeError:
                pass  # would log error here

            return a

        # if not ndarray then return as is
        return return_values
