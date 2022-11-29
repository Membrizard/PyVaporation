import typing
import json

import attr
import numpy
from pathlib import Path

from ..mixtures import Composition


# TODO: Enum or not Enum?
class CalculationType:
    polynomial: str = "polynomial"
    exponential: str = "exponential"
    logarithmic: str = "logarithmic"


@attr.s(auto_attribs=True)
class TemperatureProgram:
    coefficients: typing.List[float]
    type: str = CalculationType.polynomial

    def polynomial(self, x: float) -> float:
        """
        Calculates temperature using polynomial approximation
        :param x: parameter used for calculation of temperature
        :return: Temperature value calculated using polynomial relation defined with .coefficients
        """
        return sum(
            [self.coefficients[i] * x**i for i in range(len(self.coefficients))]
        )

    def exponential(self, x: float) -> float:
        """
        Calculates temperature using exponential approximation
        :param x: parameter used for calculation of temperature
        :return: Temperature value calculated using exponential-polynomial relation defined with .coefficients
        """
        return self.coefficients[0] * numpy.exp(
            sum(
                [
                    self.coefficients[i] * x ** (i - 1)
                    for i in range(1, len(self.coefficients))
                ]
            )
        )

    def logarithmic(self, x: float) -> float:
        """
        Calculates temperature using logarithmic approximation
        :param x: parameter used for calculation of temperature
        :return: Temperature value calculated using logarithmic-polynomial relation defined with .coefficients
        """
        return self.coefficients[0] * numpy.log(
            sum(
                [
                    self.coefficients[i] * x ** (i - 1)
                    for i in range(1, len(self.coefficients))
                ]
            )
        )

    def program(self, time):
        """
        Calculation of the Temperature based on a given temperature program
        :param time - time in hours
        :return - Temperature in K
        """
        return getattr(self, self.type)(time)


@attr.s(auto_attribs=True)
class Conditions:
    """
    A class for specification initial conditions for modelling of the test_pervaporation processes
    """

    membrane_area: float
    initial_feed_temperature: float
    initial_feed_amount: float
    initial_feed_composition: Composition
    permeate_temperature: typing.Optional[float] = None
    permeate_pressure: typing.Optional[float] = None
    temperature_program: typing.Optional[TemperatureProgram] = None

    @classmethod
    def safe_load(cls, path: typing.Union[str, Path]) -> "Conditions":
        """
        :param path: Path to a json object
        :return: Returns a Conditions object (Temperature program is ignored) from a json file
        """
        with open(path, "r") as openfile:
            # Reading from json file
            json_object = json.load(openfile)
        return Conditions(
            membrane_area=json_object["membrane_area"],
            initial_feed_temperature=json_object["initial_feed_temperature"],
            initial_feed_amount=json_object["initial_feed_amount"],
            initial_feed_composition=Composition(
                p=json_object["initial_feed_composition_value"],
                type=json_object["initial_feed_composition_type"],
            ),
            permeate_temperature=json_object["permeate_temperature"],
            permeate_pressure=json_object["permeate_pressure"],
        )

    def safe_save(self, path: typing.Union[str, Path]):
        """
        :param path: Path to a json object
        :return: Saves a Conditions object (Temperature program is ignored) to a json file
        """
        json_dict = {
            "membrane_area": self.membrane_area,
            "initial_feed_temperature": self.initial_feed_temperature,
            "initial_feed_amount": self.initial_feed_amount,
            "initial_feed_composition_value": self.initial_feed_composition.p,
            "initial_feed_composition_type": self.initial_feed_composition.type,
            "permeate_temperature": self.permeate_temperature,
            "permeate_pressure": self.permeate_pressure,
        }
        with open(path, "w") as outfile:
            json.dump(json_dict, outfile)
