import typing

import attr

import numpy
from enum import Enum
from mixture import Composition, CompositionType


@attr.s(auto_attribs=True)
class CalculationType(Enum):
    polynomial: str = 'polynomial'
    exponential: str = 'exponential'
    logarithmic: str = 'logarithmic'


@attr.s(auto_attribs=True)
class TemperatureProgram:
    coefficients: typing.List[float]
    type: CalculationType

    def program(self, time):

        def polynomial(x):
            return sum([self.coefficients[i]*x**i for i in range(len(self.coefficients))])

        def exponential(x):
            return self.coefficients[0]*numpy.exp(sum([self.coefficients[i]*x**i for i in range(len(self.coefficients))]))

        def logarithmic(x):
            return self.coefficients[0] * numpy.log(
                sum([self.coefficients[i] * x ** i for i in range(len(self.coefficients))]))
        d = {'polynomial': polynomial(time), 'exponential': exponential(time), 'logarithmic': logarithmic(time)}

        return d[self.type]


@attr.s(auto_attribs=True)
class Conditions:
    membrane_area: float
    initial_feed_temperature: float
    initial_feed_amount: float
    initial_feed_composition: Composition
    permeate_temperature: float
    temperature_program: typing.Optional[TemperatureProgram] = None



