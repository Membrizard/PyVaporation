import attr
import typing

from ..mixture import Mixture


@attr.s(auto_attribs=True)
class NonIdealExperiment:
    """
    A setting of non-ideal experiment
    """
    temperature: float  # K
    # Permeance in kg*mcm/(m2*h*kPa)
    overall_flux: typing.List[float]
    component_fluxes: typing.List[typing.List[float]]
    compositions: typing.List[float]
    mixture: Mixture
