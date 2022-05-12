import pathlib
import typing

import attr
import numpy
import yaml

from ..component import AllComponents, Component
from ..utils import Composition, NRTLParameters, R


@attr.s(auto_attribs=True)
class Mixture:
    components: typing.List[Component]
    nrtl_params: NRTLParameters

    @classmethod
    def from_dict(cls, d: typing.Mapping, all_components: AllComponents) -> "Mixture":
        return Mixture(
            components=[
                all_components.__getattribute__(component_name)
                for component_name in d["components"]
            ],
            nrtl_params=NRTLParameters(**d["nrtl_params"]),
        )

    def get_nrtl_partial_pressures(self, temperature: float, composition: Composition):
        tau = numpy.array(
            [
                self.nrtl_params.g12 / (R * temperature),
                self.nrtl_params.g21 / (R * temperature),
            ]
        )

        g_exp = numpy.exp(-tau * self.nrtl_params.alpha)

        activity_coefficients = [
            numpy.exp(
                (composition[1] ** 2)
                * (
                    tau[1]
                    * (g_exp[1] / (composition[0] + composition[1] * g_exp[1])) ** 2
                    + tau[0]
                    * g_exp[0]
                    / (composition[1] + composition[0] * g_exp[0]) ** 2
                )
            ),
            numpy.exp(
                (composition[0] ** 2)
                * (
                    tau[0]
                    * (g_exp[0] / (composition[1] + composition[0] * g_exp[0])) ** 2
                    + tau[1]
                    * g_exp[1]
                    / (composition[0] + composition[1] * g_exp[1]) ** 2
                )
            ),
        ]

        return (
            self.components[0].get_antoine_pressure(temperature)
            * activity_coefficients[0]
            * composition[0],
            self.components[1].get_antoine_pressure(temperature)
            * activity_coefficients[1]
            * composition[1],
        )


@attr.s(auto_attribs=True)
class AllMixtures:
    mixtures: typing.Mapping[str, Mixture]

    @classmethod
    def load(cls, path: pathlib.Path, all_components: AllComponents) -> "AllMixtures":
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
