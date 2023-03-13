import typing

import attr
import numpy

from ..components import Component
from ..utils import NRTLParameters, R, UNIQUACParameters
from ..activity_coefficient_models import (
    uniquac_activity_coefficient_equation,
    nrtl_activity_coefficient_equation,
    ActivityCoefficientModel,
)


def _is_in_0_to_1_range(instance: typing.Any, attribute, value: float) -> None:
    if not 0 <= value <= 1:
        raise ValueError("Give %s value is not in [0, 1] range" % value)


class CompositionType:
    """
    A class to describe type of the composition
    """

    molar: str = "molar"
    weight: str = "weight"


# class Mixture:
#     """
#     A class to represent mixtures
#     """
#
#     name: str
#     components: typing.List[Component]
#     nrtl_params: typing.Optional[NRTLParameters] = None
#     uniquac_params: typing.Optional[UNIQUACParameters] = None
#
#     def __init__(
#             self,
#             name: str,
#             components: typing.Optional[typing.List[Component]] = None,
#             nrtl_params: typing.Optional[NRTLParameters] = None,
#             uniquac_params: typing.Optional[UNIQUACParameters] = None,
#             first_component: typing.Optional[Component] = None,
#             second_component: typing.Optional[Component] = None,
#             third_component: typing.Optional[Component] = None,
#             fourth_component: typing.Optional[Component] = None,
#             fith_component: typing.Optional[Component] = None,
#     ):
#         self.name = name
#
#     def __attrs_post_init__(self):
#         if self.nrtl_params is None and self.uniquac_params is None:
#             raise ValueError(
#                 "Component Interaction parameters are required to create a mixture!"
#             )
#
#     def first_component(self):
#         return self.components[0]
#
#     def second_component(self):
#         return self.components[1]


@attr.s(auto_attribs=True)
class Mixture:
    """
    A class to represent mixtures
    """

    name: str
    first_component: typing.Optional[Component] = None
    second_component: typing.Optional[Component] = None
    components: typing.Optional[typing.List[Component]] = []
    nrtl_params: typing.Optional[NRTLParameters] = None
    uniquac_params: typing.Optional[UNIQUACParameters] = None

    def __attrs_post_init__(self):
        if not self.components:
            if self.first_component is None or self.second_component is None:
                raise ValueError(
                    "At least two Components are required to create a mixture!"
                )
            self.components.append(self.first_component)
            self.components.append(self.second_component)

        if not self.components and (self.first_component is not None or self.second_component is not None):
            raise ValueError(
                "Mixture should be created either by specifying a List of Components\
                 or by specifying first and second Component individually!"
            )

        if len(self.components) == 1 or (self.first_component is None or self.second_component is None):
            raise ValueError(
                "At least two Components are required to create a mixture!"
            )

        if self.nrtl_params is None and self.uniquac_params is None:
            raise ValueError(
                "Component Interaction parameters are required to create a mixture!"
            )

    # def first_component(self):
    #     return self.components[0]
    #
    # def second_component(self):
    #     return self.components[1]


@attr.s(auto_attribs=True)
class Composition:
    """
    A class to represent composition of the mixtures
    """

    p: float = attr.ib(validator=_is_in_0_to_1_range)
    type: str

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


def get_partial_pressures(
    temperature: float,
    mixture: Mixture,
    composition: Composition,
    calculation_type: str = ActivityCoefficientModel.NRTL,
) -> typing.Tuple[float, float]:
    """
    Calculation of partial pressures of both test_components
    :params
    temperature: temperature in K
    mixture: a mixture for which the calculation should be conducted
    composition: specified composition in mol or weight %
    :return: Partial pressures as a tuple, test_components wise in kPa
    """
    if composition.type == CompositionType.weight:
        composition = composition.to_molar(mixture=mixture)

    activity_coefficients = calculate_activity_coefficients(
        temperature=temperature,
        mixture=mixture,
        composition=composition,
        calculation_type=calculation_type,
    )
    return (
        mixture.first_component.get_vapor_pressure(temperature)
        * activity_coefficients[0]
        * composition.first,
        mixture.second_component.get_vapor_pressure(temperature)
        * activity_coefficients[1]
        * composition.second,
    )


def calculate_activity_coefficients(
    temperature: float,
    mixture: Mixture,
    composition: Composition,
    calculation_type: str = ActivityCoefficientModel.NRTL,
) -> typing.Tuple:
    """
    Calculation of activity coefficients of both test_components
    :params
    temperature: temperature in K
    mixture: a mixture for which the calculation should be conducted
    composition: specified composition in mol or weight %
    :return: activity coefficients as a tuple
    """
    if composition.type == CompositionType.weight:
        composition = composition.to_molar(mixture=mixture)

    if calculation_type == ActivityCoefficientModel.NRTL:

        if mixture.nrtl_params is None:
            raise ValueError(
                "NRTL Parameters must be specified for this type of calculation"
            )

        activity_coefficients = nrtl_activity_coefficient_equation(
            nrtl_params=mixture.nrtl_params,
            temperature=temperature,
            x=[composition.first, composition.second],
        )

        return activity_coefficients

    elif calculation_type == ActivityCoefficientModel.UNIQUAC:
        # The implementation is based on https://doi.org/10.1021/i260068a028
        if composition.first == 0:
            composition = Composition(p=0.00001, type="molar")
        if composition.second == 0:
            composition = Composition(p=0.99999, type="molar")

        if mixture.uniquac_params is None:
            raise ValueError(
                "UNIQUAC Parameters must be specified for this type of calculation"
            )
        if (
            mixture.first_component.uniquac_constants is None
            or mixture.second_component.uniquac_constants is None
        ):
            raise ValueError(
                "UNIQUAC Constants for all Components must be specified for this type of calculation"
            )

        first_component_const = mixture.first_component.uniquac_constants
        second_component_const = mixture.second_component.uniquac_constants

        binary_interaction_params = mixture.uniquac_params.binary_parameters_matrix
        # tau[i][i] = 1 !!!!
        tau = [[1, 0], [0, 1]]

        for i in range(len(binary_interaction_params)):
            for j in range(len(binary_interaction_params)):
                if binary_interaction_params[i][j] != 0:
                    tau[i][j] = numpy.exp(
                        -(
                            binary_interaction_params[i][j][0]
                            + binary_interaction_params[i][j][1] / temperature
                        )
                        / temperature
                    )

        activity_coefficients = uniquac_activity_coefficient_equation(
            r=[first_component_const.r, second_component_const.r],
            q_geometric=[
                first_component_const.q_geometric,
                second_component_const.q_geometric,
            ],
            q_interaction=[
                first_component_const.q_interaction,
                second_component_const.q_interaction,
            ],
            x=[composition.first, composition.second],
            tau=tau,
            z=mixture.uniquac_params.z,
        )

        return activity_coefficients
