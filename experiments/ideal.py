import attr
import typing

from ..component import Component


@attr.s(auto_attribs=True)
class IdealExperiment:
    """
    A setting of ideal experiment.
    """
    temperature: float  # K
    permeance: float  # kg * mcm / (m2 * h * kPa)
    component: Component
    activation_energy: typing.Optional[float] = None
