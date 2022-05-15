import typing
from pathlib import Path

import attr

from mixture import Composition, Mixture


@attr.s(auto_attribs=True)
class DiffusionCurve:
    mixture: Mixture
    membrane_name: str
    temperature_feed: float
    temperature_permeate: typing.Optional[float]
    compositions: typing.List[Composition]
    partial_flux_1: typing.List[float]
    partial_flux_2: typing.List[float]
    permeance_1: typing.Optional[typing.List[float]]
    permeance_2: typing.Optional[typing.List[float]]
    comments: typing.Optional[str]

    @property
    def separation_factor(self) -> typing.List[float]:
        return [0]

    @property
    def psi(self) -> typing.List[float]:
        return [0]

    @property
    def permeate_composition(self) -> typing.List[Composition]:
        return []

    @classmethod
    def from_csv(cls, path: typing.Union[str, Path]) -> "DiffusionCurve":
        pass
