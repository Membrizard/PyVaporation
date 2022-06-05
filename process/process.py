import typing
from pathlib import Path

import attr
import numpy

from conditions import Conditions
from mixture import Composition, Mixture
from permeance import Permeance


@attr.s(auto_attribs=True)
class ProcessModel:
    mixture: Mixture
    membrane_name: str
    isothermal: bool
    feed_temperature: typing.List[float]
    feed_composition: typing.List[Composition]
    permeate_composition: typing.List[Composition]
    permeate_temperature: typing.List[float]
    feed_mass: typing.List[float]
    partial_fluxes: typing.List[typing.Tuple[float, float]]
    permeances: typing.List[typing.Tuple[Permeance, Permeance]]
    time: typing.List[float]
    feed_evaporation_heat: typing.List[float]
    permeate_condensation_heat: typing.List[float]
    initial_conditions: Conditions
    comments: typing.Optional[str] = None

    @property
    def get_separation_factor(self) -> typing.List[float]:
        """
        :return: List of separation Factors
        """
        feed = self.feed_composition
        permeate = self.permeate_composition
        return [
            ((1 - feed[i].second) / feed[i].second)
            / ((1 - permeate[i].second) / permeate[i].second)
            for i in range(len(feed))
        ]

    @property
    def get_psi(self) -> typing.List[float]:
        """
        :return: List of Pervaporation Separation Index (PSI) values
        """
        separation_factor = self.get_separation_factor
        total_flux = [
            sum(self.partial_fluxes[i]) for i in range(len(self.partial_fluxes))
        ]
        return [numpy.multiply(total_flux, numpy.substract(separation_factor, 1))]

    @property
    def get_selectivity(self) -> typing.List[float]:
        """
        :return: List of Calculated selectivities
        """
        permeance = self.permeances
        return [permeance[i][0] / permeance[i][1] for i in range(len(permeance))]

    @classmethod
    def from_csv(cls, path: typing.Union[str, Path]) -> "ProcessModel":
        pass
