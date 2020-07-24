#!/usr/bin/env python
# coding: utf-8


class AttrDict(object):
    """Class with Dictionary-like Setting and Getting"""
    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __delitem__(self, key):
        del self.__dict__[key]

    def __contains__(self, key):
        return key in self.__dict__

    def __len__(self):
        return len(self.__dict__)

    def __repr__(self):
        return repr(self.__dict__)

    def __eq__(self, other):
        """This ensure equality between two dictionary"""
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __ne__(self, other):
        """Not necessary for Python 3"""
        return not self.__eq__(other)

    def update(self, init=None, **kwargs):
        """Update existing parameters

        Args:
            init: If init is present and has a .keys() method,
                then does: for k in init: D[k] = E[k].
                If init is present and lacks a .keys() method,
                then does: for k, v in init: D[k] = v
        """
        if isinstance(init, dict):
            self.__dict__.update(init)
        elif isinstance(init, AttrDict):
            self.__dict__.update(init.__dict__)
        else:
            self.__dict__.update(**kwargs)


class Parameter(AttrDict):
    """Parent Parameter Class"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
