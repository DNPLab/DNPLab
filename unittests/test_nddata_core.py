import pytest
import operator
import unittest
from numpy.testing import assert_array_equal
from dnplab.core.base import ABCData
from dnplab.core.coord import Coords
import numpy as np
import random

import logging

logger = logging.getLogger(__name__)

test_dims = ["x", "y", "z", "p", "q", "r"]
num_random_tests = 10


def get_random_ABCData_list(seed_axis=0, seed_data=0):
    """
    Pre-generate a list of randomized tuple (ABCData, np.array) for all the
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
    ABCData_list = []
    for random_axis in random_axes:
        dims = [axis[0] for axis in random_axis]
        coords = [axis[1] for axis in random_axis]
        shape = [coord.size for coord in coords]
        values = np.random.randn(*shape)
        ABCData_list.append((ABCData(values, dims, coords), values))
    # ABCData_list.append((ABCData.ABCData_core(), np.array([])))  # UserWarning: Github #37
    return ABCData_list


random_ABCData_list = get_random_ABCData_list(seed_axis=0, seed_data=0)
random_ABCData_list_2 = get_random_ABCData_list(seed_axis=0, seed_data=1)


@pytest.mark.filterwarnings("ignore:divide by zero")
# See https://docs.python.org/3/library/operator.html#mapping-operators-to-functions
@pytest.mark.parametrize(
    "operator", [operator.add, operator.sub, operator.mul, operator.truediv]
)
@pytest.mark.parametrize("ABCData_value_tuple", random_ABCData_list)
@pytest.mark.parametrize("number", [-1.1j, -1.1, -1, 0, 1, 1.1, 1.1j])
def test_ABCData_core_math_operators_numeric(operator, ABCData_value_tuple, number):
    ABCData, values = ABCData_value_tuple
    assert_array_equal(operator(ABCData, number).values, operator(values, number))
    assert_array_equal(operator(number, ABCData).values, operator(number, values))


@pytest.mark.parametrize(
    "operator", [operator.add, operator.sub, operator.mul, operator.truediv]
)
@pytest.mark.parametrize("i_data", range(0, len(random_ABCData_list)))
def test_ABCData_core_math_operators(operator, i_data):
    ABCData, values = random_ABCData_list[i_data]
    ABCData_2, values_2 = random_ABCData_list_2[i_data]
    ABCData_3, values_3 = operator(ABCData, ABCData_2), operator(values, values_2)
    # assert self consistent
    assert ABCData_3._self_consistent()
    # assert values equal
    assert_array_equal(ABCData_3.values, values_3)


def test_ABCData_core_math_div_by_zero():
    with pytest.warns(RuntimeWarning):
        ABCData(
            values=np.array([1, 2, 3]), coords=[np.array([3, 2, 1])], dims=["x"]
        ) / 0


class dnplab_ABCData_core_tester(unittest.TestCase):
    def setUp(self):
        self.dims = test_dims
        random.sample(test_dims, random.randint(1, len(test_dims)))

    def construct_random_data(self):
        random_dims = random.sample(test_dims, random.randint(1, len(test_dims)))

        random_coords = [np.r_[0 : random.randint(1, 6)] for dim in random_dims]
        shape = [coord.size for coord in random_coords]

        random_values = np.random.randn(*shape)
        data = ABCData(random_values, random_dims, random_coords)
        return data, random_dims, random_values, random_coords

    def test_ABCData_core_init(self):
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
        self.assertEqual(type(data.coords), Coords)

    def test_ndim(self):
        values = np.r_[1:10].reshape(3, 3)
        x = np.r_[0:3]
        y = np.r_[0:3]

        data = ABCData(values, ["x", "y"], [x, y])

        self.assertEqual(data.ndim, 2)

    def test_000_checkDimdeepcopy(self):
        # test for issue 295
        # this test doesn't really test anything but it checks that the functionality is as is, so this test fails if an error is raised
        # uses dnplab
        import dnplab as dnp

        dims = ["Average", "t2"]
        coords1 = [np.arange(0, 100), np.arange(0, 1024)]
        coords2 = [np.arange(0, 100), np.arange(0, 1024)]
        data1 = np.random.random((100, 1024))
        data2 = np.random.random((100, 1024))

        DNPObj1 = dnp.DNPData(data1, dims, coords1)
        DNPObj2 = dnp.DNPData(data2, dims, coords2)

        DNPObj1 = dnp.fourier_transform(DNPObj1, convert_to_ppm=False)
        # this should not raise an error
        try:
            DNPObj2 = dnp.fourier_transform(DNPObj2, convert_to_ppm=False)
        except ValueError as e:
            self.fail("Deepcopy of dims most likely failed! Error {0}".format(e))

    def test_001_checkUnfoldDim(self):
        # test for fold/unfold functionality, and make sure that unfold is 2d
        import dnplab as dnp

        dims = ["Average", "t2"]
        coords1 = [np.arange(0, 100), np.arange(0, 1024)]
        coords2 = [np.arange(0, 100), np.arange(0, 1024)]
        data1 = np.random.random((100, 1024))
        data2 = np.random.random(100)

        DNPObj1 = dnp.DNPData(data1, dims, coords1)
        DNPObj2 = dnp.DNPData(data2, ["Average"], [coords2[0]])

        # even only 1D objects become atleast 2d
        DNPObj2.unfold("Average")
        self.assertTrue(len(DNPObj2.values.shape) == 2)
        self.assertEqual(DNPObj2.values.shape, (100, 1))

        DNPObj1.unfold("Average")
        self.assertTrue(len(DNPObj1.values.shape) == 2)
        self.assertEqual(DNPObj1.values.shape, (100, 1024))


class dnplab_ABCData_coord_tester(unittest.TestCase):
    def setUp(self):
        self.coord_inst_a = np.r_[0:10:1]
        self.coord_inst_b = np.r_[0:1:0.005]
        self.numpy_inst = np.r_[1:2:0.25]
        self.collection_inst = Coords(
            ["a", "b", "c"], [self.coord_inst_a, self.coord_inst_b, self.numpy_inst]
        )

    def test_get_str_uses_dim(self):
        assert_array_equal(self.collection_inst.dims, ["a", "b", "c"])

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


class ABCData_numpy_implementation_test(unittest.TestCase):
    # c&p from above for simplicity
    def construct_random_data(self):
        random_dims = random.sample(test_dims, random.randint(1, len(test_dims)))

        random_coords = [np.r_[0 : random.randint(1, 6)] for dim in random_dims]
        shape = [coord.size for coord in random_coords]

        random_values = np.random.randn(*shape)
        # logger.debug('ABCData has random_dims: {0} and radom_values shape {1}'.format(random_dims,random_values))
        data = ABCData(random_values, random_dims, random_coords)
        return data, random_dims, random_values, random_coords

    def setUp(self):
        self.size = 10
        self.val = 5
        self.a = np.zeros(self.size)
        self.b = np.ones(self.size)
        self.c = self.val * np.ones(self.size)
        self.d = np.random.random((self.size, self.size))

        self.ax = np.arange(self.size)

        self.Adata = ABCData(self.a, ["dim0"], [np.arange(self.a.size)])
        self.Bdata = ABCData(self.b, ["dim0"], [np.arange(self.b.size)])
        self.Cdata = ABCData(self.c, ["dim0"], [np.arange(self.a.size)])
        self.Ddata = ABCData(self.d, dims=["1", "2"], coords=[self.ax, self.ax])

    def test_000_replaceAttr(self):
        from dnplab.core.base import _replaceClassWithAttribute

        test_tuple = (self.Adata, 1, self.a)
        test_dict = {"1": self.Bdata, "2": 5, "3": self.c}
        nargs, nkwargs = _replaceClassWithAttribute(self.Adata, test_tuple, test_dict)
        self.assertEqual(type(nargs[0]), type(self.a))

    def test_001_npSum(self):
        self.assertEqual(np.sum(self.Adata), 0)
        self.assertEqual(np.sum(self.Bdata), 10)
        self.assertEqual(np.sum(self.Cdata), 50)

    def test_002_npSin(self):
        f = np.sin
        self.assertEqual(np.sum(f(self.Adata) - f(self.a)), 0)
        self.assertEqual(np.sum(f(self.Bdata) - f(self.b)), 0)
        self.assertEqual(np.sum(f(self.Cdata) - f(self.c)), 0)
        self.assertEqual(type(f(self.Adata)), type(self.Adata))

        # check that dimensions are not deleted!
        self.assertEqual(f(self.Ddata).dims, self.Ddata.dims)
        self.assertEqual(f(self.Bdata).dims, self.Bdata.dims)

    def test_003_axisKeywordTest_nonufunc(self):
        # testing with function sum, which is often used
        b = np.sum(self.Adata)
        self.assertEqual(type(b), np.float64)

        rdata, dims, values, coords = self.construct_random_data()
        test_sum = np.sum(values, axis=0)  # corresponds to dim[0]

        logger.debug(
            "Logging rdata.shape: {0} and dims: {1}".format(
                rdata._values.shape, rdata.dims
            )
        )

        b = np.sum(rdata, axis=dims[0])

        logger.debug(
            "Logging b.shape: {0} and test_sum.shape: {1}".format(
                b._values.shape, test_sum.shape
            )
        )
        logger.debug(
            "Logging b dims: {0} and rdata.dims: {1}".format(b.dims, rdata.dims)
        )

        self.assertEqual(type(b), type(rdata))
        self.assertEqual(test_sum.shape, b.shape)
        self.assertTrue(
            np.all(np.isclose(test_sum - b._values, 0))
        )  # check for value equality
        self.assertFalse(dims[0] in b.dims)
        self.assertRaises(ValueError, np.sum, b, axis="doesnotexist")

        d=np.sum(self.Ddata,axis='2')
        self.assertEqual(d.dims,['1'])
        self.assertTrue( np.all(np.isclose(d.values-self.Ddata.values.sum(axis=-1),0)) )
