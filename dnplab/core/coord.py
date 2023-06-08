from __future__ import division
import numpy as _np
import operator
import functools
from collections import OrderedDict
from copy import deepcopy


class Coords(object):
    def __init__(self, dims, coords):
        """Object for storing axes information

        Attributes:
            dims (list): List of strings labeling dimensions
            coords (list): List of arrays for axes
        """
        # super(Coords, self).__init__()

        if not isinstance(dims, list):
            raise ValueError("dims must be a list")
        if not isinstance(coords, list):
            raise ValueError("coords must be a list")

        for index, coord in enumerate(coords):
            if isinstance(coord, (range, list)):
                coords[index] = _np.array(coord)

        if self._check_dims(dims):
            self._dims = deepcopy(dims)
        else:
            raise TypeError(
                "dims must be list of unique strings (type str), you provided types {0}".format(
                    [type(k) for k in dims]
                )
            )

        if self._check_coords(coords):
            self._coords = deepcopy(coords)
        else:
            raise TypeError(
                "coords must be list of 1d numpy arrays, you provided types {0}".format(
                    [type(k) for k in coords]
                )
            )

    def _check_dims(self, dims):
        """Verify dims is a list of str

        Args:
            dims: Object to test for valid dims

        Returns:
            bool: True if valid dims, False otherwise
        """
        # Check if dims is list
        if not isinstance(dims, list):
            return False

        # Check if all members are str
        for dim in dims:
            if not isinstance(dim, str):
                return False

        # Check for duplicates
        if len(dims) != len(set(dims)):
            return False

        return True

    def _check_coords(self, coords):
        """Check if valid coords

        Args:
            coords: Object to test

        Returns:
            bool: True if list of 1d numpy arrays, False otherwise
        """

        # Verify coords is a list
        if not isinstance(coords, list):
            return False

        # Verify each member is 1d numpy array (and not empty)
        for coord in coords:
            if not isinstance(coord, _np.ndarray):
                return False
            if (not (len(coord.shape) == 1)) or (coord.size == 0):
                return False

        return True

    def _self_consistent(self):
        return NotImplemented

    def index(self, dim):
        return self._dims.index(dim)

    def __getitem__(self, dim):
        if isinstance(dim, str):
            return self.coords[self.index(dim)]
        elif isinstance(dim, int):
            return self.coords[dim]
        else:
            raise TypeError("dim must be type str or int not: %s" % str(type(dim)))

    def __setitem__(self, dim, coord):
        if not isinstance(dim, str):
            raise TypeError("dim must be type str not %s" % str(type(dim)))
        if not isinstance(coord, (_np.ndarray)):
            raise TypeError(
                "argument must be type nddata_coord or numpy ndarray not %s"
                % str(type(coord))
            )

        # if dim already in dims, overwrite
        if dim in self.dims:
            index = self.index(dim)
            self._coords[index] = coord
        else:
            self._dims.append(dim)
            self._coords.append(coord)

    def __delitem__(self, dim):
        index = self.index(dim)
        del self._dims[index]
        del self._coords[index]

    @property
    def dims(self):
        return self._dims

    @dims.setter
    def dims(self, dims):
        if self._check_dims(dims):
            self._dims = dims
        else:
            raise TypeError("Invalid dims. Cannot set dims to {}".format(dims))

    @property
    def coords(self):
        return self._coords

    @coords.setter
    def coords(self, coords):
        if self._check_coords(coords):
            self._coords = coords
        else:
            raise TypeError("Invalid coords. Cannot set coords to {}".format(coords))

    def __repr__(self):
        return "Coords({})".format(self.coords)

    def __str__(self):
        return "dims:\n{}\ncoords:\n{}".format(self.dims, self.coords)

    def __iter__(self):
        return iter(self.coords)

    def __len__(self):
        return len(self.coords)

    @property
    def shape(self):
        return tuple(k.size for k in self.coords)

    @property
    def size(self):
        return functools.reduce(operator.mul, [len(k) for k in self.coords], 1)

    def reorder(self, dims):
        if not self._check_dims(dims):
            raise TypeError("New dims must be list of str with no duplicates")
        for dim in dims:
            if dim not in self.dims:
                raise ValueError("no such dimension: %s" % dim)

        # Add original dims to end, remove duplicates
        dims = list(OrderedDict.fromkeys(dims + self.dims))

        # New indices for dims
        new_order = [dims.index(dim) for dim in self.dims]
        permutation_order = [self.dims.index(dim) for dim in dims]

        self.dims = dims

        #        self.coords = [self.coords[x] for x in new_order]
        self.coords = [self.coords[x] for x in permutation_order]

    def pop(self, dim):
        if type(dim) == str:
            index = self.index(dim)
        else:
            index = dim
        out = self._coords.pop(index)
        self.dims.pop(index)
        return out

    def copy(self):
        return deepcopy(self)

    def reorder_index(self, new_order):
        """Reorder based on index"""
        self._coords = [self._coords[x] for x in new_order]
        self._dims = [self._dims[x] for x in new_order]

    def rename(self, dim, new_dim):
        self.dims[self.index(dim)] = new_dim

    def append(self, dim, coord):
        """Append to coords"""
        if not isinstance(dim, str):
            raise TypeError("dim must be type str not %s" % str(type(dim)))

        if not isinstance(coord, (complex, float, int, _np.ndarray)):
            raise TypeError("coord must be type numpy not %s" % str(type(coord)))

        if isinstance(coord, (list, range, float, int, complex)):
            coord = _np.array(coord)

        if dim not in self.dims:
            self._dims.append(dim)
            self._coords.append(coord)
        else:
            raise ValueError("dim already in dims, cannot append to coords")
