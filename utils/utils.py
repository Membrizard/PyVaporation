import typing

import attr
import numpy

from enum import Enum

from ..mixture import Mixture

R = 8.314462


class CompositionType(Enum):
    molar: str = "molar"
    weight: str = "weight"


@attr.s(auto_attribs=True)
class Composition:
    p: float = attr.ib(validator=lambda value: 0 <= value <= 1)  # type: ignore
    type: CompositionType

    @property
    def first(self) -> float:
        return self.p

    @property
    def second(self) -> float:
        return 1 - self.p

    def to_molar(self, mixture: Mixture) -> "Composition":
        if self.type == CompositionType.molar:
            return self
        else:
            p = (self.p / mixture.first_component.molecular_weight) / (
                    self.p / mixture.first_component.molecular_weight
                    + (1 - self.p) / mixture.second_component.molecular_weight
            )
            return Composition(p=p, type=CompositionType('weight'))

    def to_weight(self, mixture: Mixture) -> "Composition":
        if self.type == CompositionType.weight:
            return self
        else:
            p = (mixture.first_component.molecular_weight * self.p) / (
                mixture.first_component.molecular_weight * self.p
                + mixture.second_component.molecular_weight * (1 - self.p)
            )
            return Composition(p=p, type=CompositionType('weight'))


@attr.s(auto_attribs=True)
class AntoineConstants:
    a: float = attr.ib(converter=lambda value: float(value))  # type: ignore
    b: float = attr.ib(converter=lambda value: float(value))  # type: ignore
    c: float = attr.ib(converter=lambda value: float(value))  # type: ignore


@attr.s(auto_attribs=True)
class NRTLParameters:
    g12: float
    g21: float
    alpha: float


@attr.s(auto_attribs=True)
class Membrane:
    ideal_experiments: typing.List[IdealExperiment]
    non_ideal_experiments: typing.List[NonIdealExperiment]

    @property
    def ideal_experiments_names(self) -> typing.List[str]:
        return [ie.component.name for ie in self.ideal_experiments]

    # Get all the penetrants the membrane was tested for
    def get_known_penetrants(self) -> typing.List[Component]:
        return numpy.unique([
            ideal_experiment.component for ideal_experiment in self.ideal_experiments
        ])

    # Picking only Experiments related to the chosen component
    def get_penetrant_data(self, component) -> typing.List[IdealExperiment]:
        return list(
            filter(
                lambda value: value.component.name in self.ideal_experiments_names,
                self.ideal_experiments
            )
        )

    # Calculate an apparent activation energy of permeation
    def calculate_activation_energy(self,component) -> float:
        
        self.get_penetrant_data(component)
        # Calculation of Least-squares Linear Fit of Ln(Permeance) versus 1/T
        return 0

    def get_permeance(self, temperature, component) -> float:
        return 0

    def get_ideal_selectivity(self, temperature, component1, component2) -> float:
        return 0


@attr.s(auto_attribs=True)
class Pervaporation:
    membrane: Membrane
    mixture: Mixture
    conditions: Conditions
    ideal: bool = True

    def get_component_flux(self):
        pass

    def get_overall_flux(self):
        pass

    def get_separation_factor(self):
        pass

    def get_psi(self):
        pass

    def get_real_selectivity(self):
        pass
    
    def model_ideal_diffusion_curve(self):
        pass

    def model_ideal_process(self):
        pass

    def model_non_ideal_process(self):
        pass
