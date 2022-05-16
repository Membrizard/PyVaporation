import typing
from pathlib import Path

import attr

from mixture import Composition, Mixture


@attr.s(auto_attribs=True)
class DiffusionCurve:
    mixture: Mixture
    membrane_name: str
    feed_temperature: float
    permeate_temperature: typing.Optional[float]
    compositions: typing.List[Composition]
    partial_fluxes: typing.List[typing.Tuple[float, float]]
    permeances: typing.Optional[typing.List[typing.Tuple[float, float]]] = None
    comments: typing.Optional[str] = None

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
