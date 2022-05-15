import typing
from pathlib import Path

import attr

from mixture import Composition, Mixture


@attr.s(auto_attribs=True)
class DiffusionCurve:
    feed_temperature: float
    mixture: Mixture
    composition_range: typing.List[Composition]
    partial_fluxes: typing.List[typing.List[float]]
    total_flux: typing.List[float]
    separation_factor: typing.List[float]
    PSI: typing.List[float]
    permeate_temperature: typing.Optional[float] = 0

    @classmethod
    def from_csv(cls, path: typing.Union[str, Path]) -> "DiffusionCurve":
        pass


@attr.s(auto_attribs=True)
class DiffusionCurves:
    pass
