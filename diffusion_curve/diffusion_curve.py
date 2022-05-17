import typing
from pathlib import Path

import attr
import numpy

from mixture import Composition, Mixture, CompositionType, get_nrtl_partial_pressures


@attr.s(auto_attribs=True)
class DiffusionCurve:
    mixture: Mixture
    membrane_name: str
    feed_temperature: float
    permeate_temperature: typing.Optional[float]
    feed_compositions: typing.List[Composition]
    partial_fluxes: typing.List[typing.Tuple[float, float]]
    permeances: typing.Optional[typing.List[typing.Tuple[float, float]]] = None
    comments: typing.Optional[str] = None

    @property
    def get_permeate_composition(self) -> typing.List[Composition]:
        return [
            Composition(
                self.partial_fluxes[i][0] / (sum(self.partial_fluxes[i])),
                type=CompositionType.weight,
            )
            for i in range(len(self.feed_compositions))
        ]

    @property
    def get_separation_factor(self) -> typing.List[float]:
        permeate_composition = self.permeate_composition
        feed_composition = self.feed_compositions
        return [
            ((1 - feed_composition[i].p) / feed_composition[i].p)
            / ((1 - permeate_composition[i].p) / permeate_composition[i].p)
            for i in range(len(self.feed_compositions))
        ]

    @property
    def get_psi(self) -> typing.List[float]:
        total_flux = [
            sum(self.partial_fluxes[i]) for i in range(len(self.feed_compositions))
        ]
        separation_factor = self.separation_factor
        return numpy.multiply(total_flux, numpy.subtract(separation_factor, 1))

    @property
    def get_permeances(self) -> typing.List[typing.Tuple[float, float]]:
        if self.permeances is not None:
            return self.permeances
        else:
            permeate_compositions = self.permeate_composition
            feed_partial_pressures = [
                get_nrtl_partial_pressures(
                    self.feed_temperature, self.mixture, composition
                )
                for composition in self.feed_compositions
            ]
            if self.permeate_temperature is None:
                return [
                    numpy.divide(self.partial_fluxes[i], feed_partial_pressures[i])
                    for i in range(len(self.feed_compositions))
                ]
            else:
                permeate_partial_pressures = [
                    get_nrtl_partial_pressures(
                        self.permeate_temperature, self.mixture, composition
                    )
                    for composition in permeate_compositions
                ]
                return [
                    numpy.divide(
                        self.partial_fluxes[i],
                        numpy.substract(
                            feed_partial_pressures[i], permeate_partial_pressures[i]
                        ),
                    )
                    for i in range(len(self.feed_compositions))
                ]

    @property
    def get_selectivity(self) -> typing.List[float]:
        permeances = self.get_permeances
        return [
            permeances[i][0] / permeances[i][1]
            for i in range(len(self.feed_compositions))
        ]

    # TODO Add interpolations with GEKKO

    @classmethod
    def from_csv(cls, path: typing.Union[str, Path]) -> "DiffusionCurve":
        pass
