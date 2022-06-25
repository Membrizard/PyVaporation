import typing
from pathlib import Path

import attr
import numpy

from conditions import Conditions
from mixtures import Composition, Mixture, get_nrtl_partial_pressures
from optimizer import PervaporationFunction
from permeance import Permeance


@attr.s(auto_attribs=True)
class ProcessModel:
    mixture: Mixture
    membrane_name: str
    feed_temperature: typing.List[float]
    feed_composition: typing.List[Composition]
    permeate_composition: typing.List[Composition]
    permeate_temperature: typing.List[float]
    permeate_pressure: typing.List[float]
    feed_mass: typing.List[float]
    partial_fluxes: typing.List[typing.Tuple[float, float]]
    permeances: typing.List[typing.Tuple[Permeance, Permeance]]
    time: typing.List[float]
    feed_evaporation_heat: typing.List[float]
    permeate_condensation_heat: typing.List[typing.Optional[float]]
    initial_conditions: Conditions
    permeance_fits: typing.Optional[
        typing.Tuple[PervaporationFunction, PervaporationFunction]
    ] = None
    comments: typing.Optional[str] = None

    def __attrs_post_init__(self):

        for i in range(len(self.time)):
            if (
                self.permeate_pressure[i] is None
                and self.permeate_temperature[i] is not None
            ):
                self.permeate_pressure[i] = sum(
                    get_nrtl_partial_pressures(
                        self.permeate_temperature[i],
                        self.mixture,
                        self.permeate_composition[i],
                    )
                )

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
        return [
            permeance[i][0].value / permeance[i][1].value for i in range(len(permeance))
        ]

    @classmethod
    def from_csv(cls, path: typing.Union[str, Path]) -> "ProcessModel":
        pass
