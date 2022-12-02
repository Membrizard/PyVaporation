import typing

import attr
import numpy

from ..components import Component
from ..utils import NRTLParameters, R


class CompositionType:
    """
    A class to describe type of the composition
    """

    molar: str = "molar"
    weight: str = "weight"


@attr.s(auto_attribs=True)
class Mixture:
    """
    A class to represent mixtures
    """

    name: str
    first_component: Component
    second_component: Component
    nrtl_params: NRTLParameters


@attr.s(auto_attribs=True)
class Composition:
    """
    A class to represent composition of the mixtures
    """

    p: float
    type: str

    def __attrs_post_init__(self):
        if self.p < 0:
            self.p = 1e-10
        if self.p > 1:
            self.p = 1-1e-10

    @property
    def first(self) -> float:
        """
        Returns fraction of the first test_components
        """
        return self.p

    @property
    def second(self) -> float:
        """
        Returns fraction of the second test_components
        """
        return 1 - self.p

    def to_molar(self, mixture: Mixture) -> "Composition":
        """
        Converts Composition to molar %
        """
        if self.type == CompositionType.molar:
            return self
        else:
            p = (self.p / mixture.first_component.molecular_weight) / (
                self.p / mixture.first_component.molecular_weight
                + (1 - self.p) / mixture.second_component.molecular_weight
            )
            return Composition(p=p, type=CompositionType.molar)

    def to_weight(self, mixture: Mixture) -> "Composition":
        """
        Converts Composition to weight %
        """
        if self.type == CompositionType.weight:
            return self
        else:
            p = (mixture.first_component.molecular_weight * self.p) / (
                mixture.first_component.molecular_weight * self.p
                + mixture.second_component.molecular_weight * (1 - self.p)
            )
            return Composition(p=p, type=CompositionType.weight)


def get_nrtl_partial_pressures(
    temperature: float, mixture: Mixture, composition: Composition
) -> typing.Tuple[float, float]:
    """
    Calculation of partial pressures of both test_components using NRTL model
    :params
    temperature: temperature in K
    NRTL parameters:
    gij in J/mol (may also be as aij+gij/RT)
    test_mixtures: Mixture
    composition: specified composition in mol or weight %
    :return: Partial pressures as a tuple, test_components wise in kPa
    """
    if composition.type == CompositionType.weight:
        composition = composition.to_molar(mixture=mixture)

    tau = numpy.array(
        [
            (mixture.nrtl_params.a12 + mixture.nrtl_params.g12 / (R * temperature)),
            (mixture.nrtl_params.a21 + mixture.nrtl_params.g21 / (R * temperature)),
        ]
    )
    if mixture.nrtl_params.alpha21 is None:
        alphas = mixture.nrtl_params.alpha12
    else:
        alphas = [mixture.nrtl_params.alpha12, mixture.nrtl_params.alpha21]

    g_exp = numpy.exp(numpy.multiply(-tau, alphas))

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
        mixture.first_component.get_vapor_pressure(temperature)
        * activity_coefficients[0]
        * composition.first,
        mixture.second_component.get_vapor_pressure(temperature)
        * activity_coefficients[1]
        * composition.second,
    )
