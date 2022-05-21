import typing
from enum import Enum
from pathlib import Path

import attr
import numpy
import yaml

from component import AllComponents, Component
from utils import NRTLParameters, R


def _is_in_0_to_1_range(
    instance: typing.Any, attribute, value: float
) -> None:  # TODO: typing
    if not 0 <= value <= 1:
        raise ValueError("Give %s value is not in [0, 1] range" % value)


class CompositionType(Enum):
    molar: str = "molar"
    weight: str = "weight"


@attr.s(auto_attribs=True)
class Mixture:
    first_component: Component
    second_component: Component
    nrtl_params: NRTLParameters

    @classmethod
    def from_dict(cls, d: typing.Mapping, all_components: AllComponents) -> "Mixture":
        return Mixture(
            first_component=all_components.__getattribute__(d["components"][0]),
            second_component=all_components.__getattribute__(d["components"][1]),
            nrtl_params=NRTLParameters(**d["nrtl_params"]),
        )


@attr.s(auto_attribs=True)
class AllMixtures:
    mixtures: typing.Mapping[str, Mixture]

    @classmethod
    def load(
        cls, path: typing.Union[str, Path], all_components: AllComponents
    ) -> "AllMixtures":
        with open(path, "r") as handle:
            _mixtures = yaml.load(handle, Loader=yaml.FullLoader)

        output = AllMixtures(
            mixtures={
                name: Mixture.from_dict(value, all_components)
                for name, value in _mixtures.items()
            }
        )

        for name, mixture in output.mixtures.items():
            setattr(output, name, mixture)

        return output


@attr.s(auto_attribs=True)
class Composition:
    p: float = attr.ib(validator=_is_in_0_to_1_range)
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
            return Composition(p=p, type=CompositionType("molar"))

    def to_weight(self, mixture: Mixture) -> "Composition":
        if self.type == CompositionType.weight:
            return self
        else:
            p = (mixture.first_component.molecular_weight * self.p) / (
                mixture.first_component.molecular_weight * self.p
                + mixture.second_component.molecular_weight * (1 - self.p)
            )
            return Composition(p=p, type=CompositionType("weight"))


def get_nrtl_partial_pressures(
    temperature: float, mixture: Mixture, composition: Composition
) -> typing.Tuple[float, float]:

    if composition.type == CompositionType.weight:
        composition = composition.to_molar(mixture=mixture)

    tau = numpy.array(
        [
            mixture.nrtl_params.g12 / (R * temperature),
            mixture.nrtl_params.g21 / (R * temperature),
        ]
    )

    g_exp = numpy.exp(-tau * mixture.nrtl_params.alpha)

    activity_coefficients = [
        numpy.exp(
            (composition.second**2)
            * (
                tau[1]
                * (g_exp[1] / (composition.first + composition.second * g_exp[1])) ** 2
                + tau[0]
                * g_exp[0]
                / (composition.second + composition.first * g_exp[0]) ** 2
            )
        ),
        numpy.exp(
            (composition.first**2)
            * (
                tau[0]
                * (g_exp[0] / (composition.second + composition.first * g_exp[0])) ** 2
                + tau[1]
                * g_exp[1]
                / (composition.first + composition.second * g_exp[1]) ** 2
            )
        ),
    ]

    return (
        mixture.first_component.get_antoine_pressure(temperature)
        * activity_coefficients[0]
        * composition.first,
        mixture.second_component.get_antoine_pressure(temperature)
        * activity_coefficients[1]
        * composition.second,
    )
