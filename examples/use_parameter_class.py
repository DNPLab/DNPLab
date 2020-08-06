"""This example shows ways to use Parameter class"""
from dnpLab.dnpHydration import Parameter

"""1. Create a Child Parameter Class"""
class MyParam(Parameter):
    pass

"""2. Getting and Setting Parameters"""
param = MyParam()

# Setting parameter: Class way
param.egg = 1

# Setting parameter: Dictionary way
param['egg'] = 2

# Getting parameter: Class or Dictionary way
print(param.egg, param['egg'], 'should be 2 2')

# Both ways called the identical function.


"""3. Building Parameter from existing parameter"""
param1 = MyParam()
param1.egg = 1

param2 = MyParam(param1)

print(param1.egg, param2.egg, 'should be 1 1')

# Or from existing dictionary
param3 = MyParam({'egg': 1})
print(param1.egg, param2.egg, param3.egg, 'should be 1, 1, 1')

# Or from casual keywords
param4 = MyParam(egg=1, ham=2)
print(param1.egg, param2.egg, param3.egg, param4.egg, 'should be 1 1 1 1')

"""4. Change existing parameter"""
param1 = MyParam(egg=1, ham=2)

# Class way
param1.egg += 1

# Dictionary way
param1.egg += 1

print(param1.egg, 'should be 3')

# Batch update
param1.update(egg=4)
param1.update({'ham':10})
param1.update(lettus=32, tomato=64)

print(param1.egg, param1.ham, param1.lettus, param1.tomato, 'should be 4 10 32 64')

"""5. Equality"""
basic1 = MyParam(egg=1)
basic2 = MyParam(egg=1)
delux = MyParam(egg=10)

print(basic1 == basic2, 'should be True')
print(basic1 == delux, 'should be False')
