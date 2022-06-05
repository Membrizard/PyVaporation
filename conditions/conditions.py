import typing

import attr
from pandas.compat.numpy import function

from mixture import Composition


@attr.s(auto_attribs=True)
class TemperatureProgram:
    coefs: typing.List[float]

    def temperature(self, time):
        return time*self.coefs


@attr.s(auto_attribs=True)
class Conditions:
    membrane_area: float
    initial_feed_temperature: float
    initial_feed_amount: float
    initial_feed_composition: Composition
    permeate_temperature: float
    temperature_program: typing.Optional[TemperatureProgram] = None
