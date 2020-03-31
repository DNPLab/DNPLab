# coding: utf-8

from unittest import TestCase

from odnpLab.parameter import AttrDict, Parameter


class TestAttrDict(TestCase):
    def test_should_init_with_one_dict(self):
        my_dict = AttrDict({'eggs': 42, 'spam': 'ham'})
        self.assertEquals(my_dict.eggs, 42)
        self.assertEquals(my_dict['eggs'], 42)
        self.assertEquals(my_dict.spam, 'ham')
        self.assertEquals(my_dict['spam'], 'ham')

    def test_should_init_with_another_attrdict(self):
        one_dict = AttrDict({'eggs': 42, 'spam': 'ham'})
        another_dict = AttrDict(one_dict)
        self.assertEqual(one_dict, another_dict)

    def test_should_not_change_values_by_inited_dict(self):
        base = {'eggs': 42, 'spam': 'ham'}
        my_dict = AttrDict(base)
        base['eggs'] = 123
        self.assertEquals(my_dict.eggs, 42)
        self.assertEquals(my_dict['eggs'], 42)
        self.assertEquals(my_dict.spam, 'ham')
        self.assertEquals(my_dict['spam'], 'ham')

    def test_should_not_change_values_by_update(self):
        base = {'eggs': 42, 'spam': 'ham'}
        my_dict = AttrDict(base)
        new = {'cheese': 10}
        new_dict = AttrDict(new)
        my_dict.update(new_dict)
        self.assertEqual(my_dict.eggs, 42)
        self.assertEqual(my_dict.cheese, 10)

    def test_get_item(self):
        my_dict = AttrDict()
        my_dict.test = 123
        self.assertEquals(my_dict.test, 123)
        self.assertEquals(my_dict['test'], 123)

    def test_set_item(self):
        my_dict = AttrDict()
        my_dict['test'] = 123
        self.assertEquals(my_dict['test'], 123)
        self.assertEquals(my_dict.test, 123)

    def test_del_attr(self):
        my_dict = AttrDict()
        my_dict['test'] = 123
        my_dict['python'] = 42
        del my_dict['test']
        del my_dict.python
        with self.assertRaises(KeyError):
            temp = my_dict['test']
        with self.assertRaises(AttributeError):
            temp = my_dict.python

    def test_in_should_work_like_in_dict(self):
        my_dict = AttrDict()
        my_dict['test'] = 123
        self.assertIn('test', my_dict)
        self.assertNotIn('bla', my_dict)

    def test_len_should_work_like_in_dict(self):
        my_dict = AttrDict()
        my_dict['test'] = 123
        my_dict.python = 42
        self.assertEquals(len(my_dict), 2)

    def test_repr(self):
        my_dict = AttrDict()
        my_dict['test'] = 123
        my_dict.python = 42
        self.assertEquals(repr(my_dict), "{'test': 123, 'python': 42}")

    def test_equal(self):
        my_dict = AttrDict({'test': 123})
        my_dict_2 = AttrDict({'test': 123})
        self.assertEqual(my_dict, my_dict_2)


class TestParameter(TestCase):
    def test_update_using_kwargs(self):
        my_param = Parameter()
        my_param.update(test=123, test2=321)
        self.assertEqual(my_param.test, 123)
        self.assertEqual(my_param.test2, 321)

    def test_update_using_dict(self):
        my_param = Parameter()
        my_param.update({'test': 123, 'test2': 321})
        self.assertEqual(my_param.test, 123)
        self.assertEqual(my_param.test2, 321)
