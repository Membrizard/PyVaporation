import typing

import attr
import numpy

from ..components import Component
from ..utils import NRTLParameters, R, UNIQUACParameters


def _is_in_0_to_1_range(instance: typing.Any, attribute, value: float) -> None:
    if not 0 <= value <= 1:
        raise ValueError("Give %s value is not in [0, 1] range" % value)


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
    nrtl_params: typing.Optional[NRTLParameters] = None
    uniquac_params: typing.Optional[UNIQUACParameters] = None

    def __attrs_post_init__(self):
        if self.nrtl_params is None and self.uniquac_params is None:
            raise ValueError(
                "Component Interaction parameters are required to create a mixture!"
            )


class ActivityCoefficientModel:
    """
    Class to represent types of activity coefficient model used for calculation of partial pressures
    """

    NRTL: str = "NRTL"
    UNIQUAC: str = "UNIQUAC"


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

    activity_coefficients = calculate_activity_coefficients(temperature=temperature,
                                                            mixture=mixture,
                                                            composition=composition,
                                                            calculation_type=calculation_type)
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
) -> typing.Tuple[float, float]:
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
            raise ValueError("NRTL Parameters must be specified for this type of calculation")

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

        activity_coefficients = (
            numpy.exp(
                (composition.second**2)
                * (
                    tau[1]
                    * (g_exp[1] / (composition.first + composition.second * g_exp[1]))
                    ** 2
                    + tau[0]
                    * g_exp[0]
                    / (composition.second + composition.first * g_exp[0]) ** 2
                )
            ),
            numpy.exp(
                (composition.first**2)
                * (
                    tau[0]
                    * (g_exp[0] / (composition.second + composition.first * g_exp[0]))
                    ** 2
                    + tau[1]
                    * g_exp[1]
                    / (composition.first + composition.second * g_exp[1]) ** 2
                )
            ),
        )

        return activity_coefficients

    elif calculation_type == ActivityCoefficientModel.UNIQUAC:
        # The implementation is based on https://doi.org/10.1021/i260068a028
        if composition.first == 0:
            composition = Composition(p=0.00001, type="molar")
        if composition.second == 0:
            composition = Composition(p=0.99999, type="molar")

        if mixture.uniquac_params is None:
            raise ValueError("UNIQUAC Parameters must be specified for this type of calculation")
        if mixture.first_component.uniquac_constants is None or mixture.second_component.uniquac_constants is None:
            raise ValueError("UNIQUAC Constants for all Components must be specified for this type of calculation")

        first_component_const = mixture.first_component.uniquac_constants
        second_component_const = mixture.second_component.uniquac_constants

        binary_interaction_params = mixture.uniquac_params.binary_parameters_matrix
        tau = [[0, 0], [0, 0]]

        for i in range(len(binary_interaction_params)):
            for j in range(len(binary_interaction_params)):
                if binary_interaction_params[i][j] != 0:
                    tau[i][j] = (
                        numpy.exp(-(binary_interaction_params[i][j][0]
                                    + binary_interaction_params[i][j][1] / temperature)
                                  / temperature))

        gammas = activity_coefficient_equation(
            r=[first_component_const.r,
               second_component_const.r],
            q_geometric=[first_component_const.q_geometric,
                         second_component_const.q_geometric],
            q_interaction=[first_component_const.q_interaction,
                           second_component_const.q_interaction],
            x=[composition.first, composition.second],
            tau=tau,
            z=mixture.uniquac_params.z,
        )

        # phi_sum = (
        #     composition.first * mixture.first_component.uniquac_constants.r
        #     + composition.second * mixture.second_component.uniquac_constants.r
        # )
        # phi_1 = (
        #     composition.first * mixture.first_component.uniquac_constants.r / phi_sum
        # )
        # phi_2 = (
        #     composition.second * mixture.second_component.uniquac_constants.r / phi_sum
        # )
        #
        # theta_sum_geometric = (
        #     composition.first * mixture.first_component.uniquac_constants.q_geometric
        #     + composition.second
        #     * mixture.second_component.uniquac_constants.q_geometric
        # )
        #
        # theta_1_geometric = (
        #     composition.first
        #     * mixture.first_component.uniquac_constants.q_geometric
        #     / theta_sum_geometric
        # )
        # theta_2_geometric = (
        #     composition.second
        #     * mixture.second_component.uniquac_constants.q_geometric
        #     / theta_sum_geometric
        # )
        #
        # theta_sum_interaction = (
        #     composition.first * mixture.first_component.uniquac_constants.q_interaction
        #     + composition.second
        #     * mixture.second_component.uniquac_constants.q_interaction
        # )
        #
        # theta_1_interaction = (
        #     composition.first
        #     * mixture.first_component.uniquac_constants.q_interaction
        #     / theta_sum_interaction
        # )
        # theta_2_interaction = (
        #     composition.second
        #     * mixture.second_component.uniquac_constants.q_interaction
        #     / theta_sum_interaction
        # )
        #
        # l_1 = mixture.uniquac_params.z / 2 * (
        #     mixture.first_component.uniquac_constants.r
        #     - mixture.first_component.uniquac_constants.q_geometric
        # ) - (mixture.first_component.uniquac_constants.r - 1)
        #
        # l_2 = mixture.uniquac_params.z / 2 * (
        #     mixture.second_component.uniquac_constants.r
        #     - mixture.second_component.uniquac_constants.q_geometric
        # ) - (mixture.second_component.uniquac_constants.r - 1)
        #
        # a_12 = mixture.uniquac_params.alpha_12 + mixture.uniquac_params.beta_12/temperature
        # a_21 = mixture.uniquac_params.alpha_21 + mixture.uniquac_params.beta_21/temperature
        #
        # tau_12 = numpy.exp(-a_12/temperature)
        # tau_21 = numpy.exp(-a_21/temperature)
        #
        # gamma_1 = numpy.exp(
        #     numpy.log(phi_1 / composition.first)
        #     + mixture.uniquac_params.z
        #     / 2
        #     * mixture.first_component.uniquac_constants.q_geometric
        #     * numpy.log(theta_1_geometric / phi_1)
        #     + phi_2
        #     * (
        #         l_1
        #         - mixture.first_component.uniquac_constants.r
        #         / mixture.second_component.uniquac_constants.r
        #         * l_2
        #     ) - mixture.first_component.uniquac_constants.q_interaction
        #     * numpy.log(theta_1_interaction + theta_2_interaction * tau_21)
        #     + theta_2_interaction * mixture.first_component.uniquac_constants.q_interaction
        #     * (tau_21 / (theta_1_interaction+theta_2_interaction * tau_21)
        #        - tau_12 / (theta_2_interaction + theta_1_interaction * tau_12))
        # )
        #
        # gamma_2 = numpy.exp(
        #     numpy.log(phi_2 / composition.second)
        #     + mixture.uniquac_params.z
        #     / 2
        #     * mixture.second_component.uniquac_constants.q_geometric
        #     * numpy.log(theta_2_geometric / phi_2)
        #     + phi_1
        #     * (
        #             l_2
        #             - mixture.second_component.uniquac_constants.r
        #             / mixture.first_component.uniquac_constants.r
        #             * l_1
        #     ) - mixture.second_component.uniquac_constants.q_interaction
        #     * numpy.log(theta_2_interaction + theta_1_interaction * tau_12)
        #     + theta_1_interaction * mixture.second_component.uniquac_constants.q_interaction
        #     * (tau_12 / (theta_2_interaction + theta_1_interaction * tau_21)
        #        - tau_12 / (theta_1_interaction + theta_2_interaction * tau_12))
        # )

        return gammas[0], gammas[1]


def activity_coefficient_equation(r: typing.List[float],
                                  q_geometric: typing.List[float],
                                  q_interaction: typing.List[float],
                                  x: typing.List[float],
                                  tau: typing.List[typing.List[float]],
                                  z: int):
    """
        Calculation of activity coefficients for multicomponent mixture using UNIQUAC equation
        Based on https://doi.org/10.1021/i260068a028 page 559
        :params
        r: A List of all r component UNIQUAC constants
        q_geometric: A List of all q_geometric component UNIQUAC constants
        q_interaction: A List of all q_interaction component UNIQUAC constants
        x: A list of molar fractions of each component in the Mixture
        tau: A matrix of components' binary interaction parameters
        z: UNIQUAC interaction parameter
        :return: activity coefficients as a List
    """

    phi_sum = numpy.sum(numpy.multiply(r, x))
    phi = [r[i]*x[i]/phi_sum for i in range(len(x))]

    theta_geometric_sum = numpy.sum(numpy.multiply(q_geometric, x))
    theta_geometric = [q_geometric[i]*x[i]/theta_geometric_sum for i in range(len(x))]

    theta_interaction_sum = numpy.sum(numpy.multiply(q_interaction, x))
    theta_interaction = [q_interaction[i] * x[i] / theta_interaction_sum for i in range(len(x))]

    l = [z/2*(r[i]-q_geometric[i])-(r[i]-1) for i in range(len(x))]

    gamma = []

    for i in range(len(x)):

        # Calculating first term in residual contribution
        interaction_term_1 = numpy.sum(numpy.multiply(theta_interaction, numpy.transpose(tau)[i]))

        # Calculating second term in residual contribution
        interaction_term_2 = 0
        interaction_term_2_sum = 0

        for j in range(len(x)):

            for k in range(len(x)):
                interaction_term_2_sum += theta_interaction[j]*tau[k][j]

            interaction_term_2 += theta_interaction[j]*tau[i][j]/interaction_term_2_sum

        gamma.append(numpy.exp(
            numpy.log(phi[i]/x[i])
            + z/2 * q_geometric[i] * numpy.log(theta_geometric[i] / phi[i])
            + l[i]
            - phi[i] / x[i] * numpy.sum(numpy.multiply(x, l))
            - q_interaction[i] * (
                + numpy.log(interaction_term_1)
                + interaction_term_2
                - 1)
         )
        )

    return gamma

