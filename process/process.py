import typing
import numpy
from pathlib import Path

import attr

from mixture import Composition, Mixture
from conditions import Conditions


@attr.s(auto_attribs=True)
class ProcessModel:
    mixture: Mixture
    membrane_name: str
    isothermal: bool
    feed_temperature: typing.List[float]
    permeate_temperature: typing.List[float]
    feed_composition: typing.List[Composition]
    permeate_composition: typing.List[Composition]
    feed_mass: typing.List[float]
    partial_fluxes: typing.List[typing.Tuple[float, float]]
    permeances: typing.List[typing.Tuple[float, float]]
    time: typing.List[float]
    feed_evaporation_heat: typing.List[float]
    permeate_condensation_heat: typing.List[float]
    initial_condtioins: Conditions
    IsTimeDefined: True
    comments: typing.Optional[str]

    @property
    def to_dimensionless_length(self):
        return [0]

    @property
    def get_separation_factor(self) -> typing.List[float]:
        feed = self.feed_composition
        permeate = self.permeate_composition
        return [((1 - feed.p) / feed.p) / ((1 - permeate.p) / permeate.p)]

    @property
    def get_psi(self) -> typing.List[float]:
        separation_factor = self.get_separation_factor
        total_flux = [
            sum(self.partial_fluxes[i]) for i in range(len(self.partial_fluxes))
        ]
        return [numpy.multiply(total_flux, numpy.substract(separation_factor, 1))]

    @property
    def get_selectivity(self) -> typing.List[float]:
        permeance = self.permeances
        return [permeance[i][0]/permeance[i][1] for i in range(len(permeance))]

    @classmethod
    def from_csv(cls, path: typing.Union[str, Path]) -> "ProcessModel":
        pass
