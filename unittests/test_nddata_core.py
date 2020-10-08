import pytest
import operator
import unittest
from numpy.testing import assert_array_equal
from dnplab.core import nddata
import numpy as np
import random
random.seed(0)

from dnplab.core.nddata_coord import nddata_coord, nddata_coord_collection

test_dims = ["x", "y", "z", "p", "q", "r"]
num_random_tests = 10


@pytest.mark.filterwarnings("ignore:divide by zero")
# See https://docs.python.org/3/library/operator.html#mapping-operators-to-functions
@pytest.mark.parametrize("operator", [operator.add, operator.sub, operator.mul, operator.truediv])
@pytest.mark.parametrize("random_axis", [
    random.sample(list(zip(test_dims, [np.r_[0 : random.randint(1, 6)] for i in test_dims])), 3)
    for _ in range(3)
])
@pytest.mark.parametrize("number", [-1.1j, -1.1, -1, 0, 1, 1.1, 1.1j])
def test_nddata_core_math_operators_numeric(operator, random_axis, number):
    dims = [axis[0] for axis in random_axis]
    coords = [axis[1] for axis in random_axis]
    shape = [coord.size for coord in coords]
    values = np.random.randn(*shape)
    data = nddata.nddata_core(values, dims, coords)

    assert_array_equal(operator(data, number).values, operator(values, number))
    assert_array_equal(operator(number, data).values, operator(number, values))


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

    def test_nddata_core_add(self):
        for ix in range(num_random_tests):
            random_coords = [np.r_[0 : random.randint(1, 6)] for dim in test_dims]

            random_axis = list(zip(test_dims, random_coords))

            random_axis1 = random.sample(random_axis, 3)
            random_axis2 = random.sample(random_axis, 3)

            dims1 = [axis[0] for axis in random_axis1]
            coords1 = [axis[1] for axis in random_axis1]
            shape1 = [coord.size for coord in coords1]
            values1 = np.random.randn(*shape1)
            data1 = nddata.nddata_core(values1, dims1, coords1)

            dims2 = [axis[0] for axis in random_axis2]
            coords2 = [axis[1] for axis in random_axis2]
            shape2 = [coord.size for coord in coords2]
            values2 = np.random.randn(*shape2)
            data2 = nddata.nddata_core(values2, dims2, coords2)

            data = data1 + data2

            self.assertTrue(data._self_consistent())
            assert_array_equal((data1 + 1).values, values1 + 1)
            assert_array_equal((data1 + 1.0).values, values1 + 1.0)
            assert_array_equal((data1 + 1.0j).values, values1 + 1.0j)


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
