import typing

import attr

from ..component import Component
from ..conditions import Conditions
from ..membrane import Membrane
from ..mixture import Mixture

R = 8.314462


@attr.s(auto_attribs=True)
class Composition:
    # TODO: from 2 to n
    p: float = attr.ib(validator=lambda value: 0 <= value <= 1)  # type: ignore

    def __getitem__(self, item: int):
        if item == 0:
            return self.p
        elif item == 1:
            return 1 - self.p
        else:
            raise ValueError("Index %s out of range" % item)


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


# Experiments for Ideal modelling
@attr.s(auto_attribs=True)
class IdealExperiment:
    temperature: float
    # Permeance in kg*mcm/(m2*h*kPa)
    permeance: float
    component: Component
    activation_energy: typing.Optional[float] = None


# Experiments for Non-ideal modelling
@attr.s(auto_attribs=True)
class NonIdealExperiment:
    temperature: float
    # Permeance in kg*mcm/(m2*h*kPa)
    overall_flux: typing.List[float]
    component_fluxes: typing.List[typing.List[float]]
    compositions: typing.List[float]
    mixture: Mixture


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
