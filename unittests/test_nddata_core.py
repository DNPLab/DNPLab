import pytest
import operator
import unittest
from numpy.testing import assert_array_equal
from dnplab.core import nddata
import numpy as np
import random

from dnplab.core.nddata_coord import nddata_coord, nddata_coord_collection

test_dims = ["x", "y", "z", "p", "q", "r"]
num_random_tests = 10


def get_random_nddata_list(seed_axis=0, seed_data=0):
    """
    Pre-generate a list of randomized tuple (nddata, np.array) for all the
    test cases here
    """
    random.seed(seed_axis)
    random_axes = [
        random.sample(
            list(zip(test_dims, [np.r_[0 : random.randint(1, 6)] for _ in test_dims])),
            3,
        )
        for _ in range(3)
    ]
    random.seed(seed_data)
    nddata_list = []
    for random_axis in random_axes:
        dims = [axis[0] for axis in random_axis]
        coords = [axis[1] for axis in random_axis]
        shape = [coord.size for coord in coords]
        values = np.random.randn(*shape)
        nddata_list.append((nddata.nddata_core(values, dims, coords), values))
    # nddata_list.append((nddata.nddata_core(), np.array([])))  # UserWarning: Github #37
    return nddata_list


random_nddata_list = get_random_nddata_list(seed_axis=0, seed_data=0)
random_nddata_list_2 = get_random_nddata_list(seed_axis=0, seed_data=1)


@pytest.mark.filterwarnings("ignore:divide by zero")
# See https://docs.python.org/3/library/operator.html#mapping-operators-to-functions
@pytest.mark.parametrize(
    "operator", [operator.add, operator.sub, operator.mul, operator.truediv]
)
@pytest.mark.parametrize("nddata_value_tuple", random_nddata_list)
@pytest.mark.parametrize("number", [-1.1j, -1.1, -1, 0, 1, 1.1, 1.1j])
def test_nddata_core_math_operators_numeric(operator, nddata_value_tuple, number):
    nddata, values = nddata_value_tuple
    assert_array_equal(operator(nddata, number).values, operator(values, number))
    assert_array_equal(operator(number, nddata).values, operator(number, values))


@pytest.mark.parametrize(
    "operator", [operator.add, operator.sub, operator.mul, operator.truediv]
)
@pytest.mark.parametrize("i_data", range(0, len(random_nddata_list)))
def test_nddata_core_math_operators(operator, i_data):
    nddata, values = random_nddata_list[i_data]
    nddata_2, values_2 = random_nddata_list_2[i_data]
    nddata_3, values_3 = operator(nddata, nddata_2), operator(values, values_2)
    # assert self consistent
    assert nddata_3._self_consistent()
    # assert values equal
    assert_array_equal(nddata_3.values, values_3)


def test_nddata_core_math_div_by_zero():
    with pytest.warns(RuntimeWarning):
        nddata.nddata_core(
            values=np.array([1, 2, 3]), coords=[np.array([3, 2, 1])], dims=["x"]
        ) / 0


class dnplab_nddata_core_tester(unittest.TestCase):
    def setUp(self):
        self.dims = test_dims
        random.sample(test_dims, random.randint(1, len(test_dims)))

    def construct_random_data(self):
        random_dims = random.sample(test_dims, random.randint(1, len(test_dims)))

        random_coords = [np.r_[0 : random.randint(1, 6)] for dim in random_dims]
        shape = [coord.size for coord in random_coords]

        random_values = np.random.randn(*shape)
        data = nddata.nddata_core(random_values, random_dims, random_coords)
        return data, random_dims, random_values, random_coords

    def test_nddata_core_init(self):
        for ix in range(num_random_tests):
            (
                data,
                random_dims,
                random_values,
                random_coords,
            ) = self.construct_random_data()
            self.assertTrue(data._self_consistent())
            assert_array_equal(data.values, random_values)
            for ix, dim in enumerate(random_dims):
                assert_array_equal(data.coords[dim], random_coords[ix])
            self.assertListEqual(data.dims, random_dims)

    def test_coord_type(self):
        data, _, _, _ = self.construct_random_data()
        self.assertEqual(type(data.coords), nddata_coord_collection)

    def test_ndim(self):
        values = np.r_[1:10].reshape(3, 3)
        x = np.r_[0:3]
        y = np.r_[0:3]

        data = nddata.nddata_core(values, ["x", "y"], [x, y])

        self.assertEqual(data.ndim, 2)


class dnplab_nddata_coord_tester(unittest.TestCase):
    def setUp(self):
        self.coord_inst_a = nddata_coord("a", slice(0, 10, 1))
        self.coord_inst_b = nddata_coord("b", slice(0, 1, 50e-3))
        self.numpy_inst = np.r_[1:2:0.25]
        self.collection_inst = nddata_coord_collection(
            ["a", "b", "c"], [self.coord_inst_a, self.coord_inst_b, self.numpy_inst]
        )

    def test_get_str_uses_dim(self):
        assert_array_equal(self.collection_inst.dims, ["a", "b", "c"])

    def test_str_rep(self):
        self.assertEqual(str(self.collection_inst["a"]), r"'a':[0 1 2 3 4 5 6 7 8 9]")

    def test_get_int_uses_index(self):
        assert_array_equal(self.collection_inst[0], self.coord_inst_a)

    def test_like_dict(self):
        assert_array_equal(self.collection_inst["a"], self.coord_inst_a)

    def test_reorder_only_1_dim(self):
        self.collection_inst.reorder(["c"])
        assert_array_equal(self.collection_inst.dims, ["c", "a", "b"])

    def test_reorder_all_3_dims(self):
        self.collection_inst.reorder(["b", "c", "a"])
        assert_array_equal(self.collection_inst.dims, ["b", "c", "a"])

    def test_rename(self):
        self.collection_inst.rename("a", "new_a")
        assert_array_equal(self.collection_inst.dims, ["new_a", "b", "c"])

    def test_rename_child_rep_changes(self):
        self.collection_inst.rename("a", "new_a")
        self.assertEqual(
            str(self.collection_inst["new_a"]), r"'new_a':[0 1 2 3 4 5 6 7 8 9]"
        )
