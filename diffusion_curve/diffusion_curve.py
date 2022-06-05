import typing
from pathlib import Path

import attr
import numpy

from mixture import Composition, CompositionType, Mixture, get_nrtl_partial_pressures


@attr.s(auto_attribs=True)
class DiffusionCurve:
    """
    The class for Diffusion Curves - dependency of selective-transport properties
    on feed composition, and process temperature
    """

    mixture: Mixture
    membrane_name: str
    feed_temperature: float
    permeate_temperature: typing.Optional[float]
    feed_compositions: typing.List[Composition]
    partial_fluxes: typing.List[typing.Tuple[float, float]]
    permeances: typing.Optional[typing.List[typing.Tuple[float, float]]] = None
    comments: typing.Optional[str] = None

    def __attrs_post_init__(self):
        """
        Calculation of Permeance values if they are not specified
        Permeance values are in kg/(m2*h*kPa)
        If needed, Permeance values may be converted using Permeance.converter()
        :return a list of Permeances for each component tuple(Pi,Pj) at each concentration
        """
        if self.permeances is None:
            permeate_compositions = self.permeate_composition
            feed_partial_pressures = [
                get_nrtl_partial_pressures(
                    self.feed_temperature, self.mixture, composition
                )
                for composition in self.feed_compositions
            ]
            if self.permeate_temperature is None:
                self.permeances = [
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
                self.permeances = [
                    numpy.divide(
                        self.partial_fluxes[i],
                        numpy.substract(
                            feed_partial_pressures[i], permeate_partial_pressures[i]
                        ),
                    )
                    for i in range(len(self.feed_compositions))
                ]

    @property
    def permeate_composition(self) -> typing.List[Composition]:
        """
        Calculation of permeate compositions at each concentration
        :return - List of permeate weight compositions
        """
        return [
            Composition(
                self.partial_fluxes[i][0] / (sum(self.partial_fluxes[i])),
                type=CompositionType("weight"),
            )
            for i in range(len(self.feed_compositions))
        ]

    @property
    def get_separation_factor(self) -> typing.List[float]:
        """
        Calculation of separation factor at each concentration
        :return - List of separation factors
        """
        permeate_composition = self.permeate_composition
        feed_composition = self.feed_compositions
        return [
            ((1 - feed_composition[i].second) / feed_composition[i].p)
            / ((1 - permeate_composition[i].second) / permeate_composition[i].p)
            for i in range(len(self.feed_compositions))
        ]

    @property
    def get_psi(self) -> typing.List[float]:
        """
        Calculation of Pervaporation Separation Index (PSI) at each concentration
        :return - List of PSI
        """
        total_flux = [
            sum(self.partial_fluxes[i]) for i in range(len(self.feed_compositions))
        ]
        separation_factor = self.get_separation_factor
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
    def get_selectivity(
        self,
        calculation_type: typing.Optional[CompositionType] = CompositionType("molar"),
    ) -> typing.List[float]:
        """
        Calculation of selectivity at each concentration
        :param calculation_type - specifies which type is of selectivity is calculated,
         based on weight or molar permeances
        :return - List of selectivity (Permeances are in kg/(m2*h*kPa) by default)
        """
        permeances = self.permeances
        if calculation_type == CompositionType("weight"):
            return [
                permeances[i][0] / permeances[i][1]
                for i in range(len(self.feed_compositions))
            ]
        else:
            return [
                permeances[i][0]
                / permeances[i][1]
                * self.mixture.second_component.molecular_weight
                / self.mixture.first_component.molecular_weight
                for i in range(len(self.feed_compositions))
            ]

    @classmethod
    def from_csv(cls, path: typing.Union[str, Path]) -> "DiffusionCurve":
        pass


@attr.s(auto_attribs=True)
class DiffusionCurves:
    diffusion_curves: typing.List[DiffusionCurve]
