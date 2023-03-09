import numpy

class ActivityCoefficientModel:
    """
    Class to represent types of activity coefficient model used for calculation of partial pressures
    """

    NRTL: str = "NRTL"
    UNIQUAC: str = "UNIQUAC"


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
        # tau[i][i] = 1 !!!!
        tau = [[1, 0], [0, 1]]

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

        return gammas[0], gammas[1]


def uniquac_activity_coefficient_equation(r: typing.List[float],
                                  q_geometric: typing.List[float],
                                  q_interaction: typing.List[float],
                                  x: typing.List[float],
                                  tau: typing.List[typing.List[float]],
                                  z: int) -> typing.Tuple[float]:
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

    phi_sum = numpy.dot(r, x)
    phi = [r[i]*x[i]/phi_sum for i in range(len(x))]

    theta_geometric_sum = numpy.dot(q_geometric, x)
    theta_geometric = [q_geometric[i]*x[i]/theta_geometric_sum for i in range(len(x))]

    theta_interaction_sum = numpy.dot(q_interaction, x)
    theta_interaction = [q_interaction[i] * x[i] / theta_interaction_sum for i in range(len(x))]

    l = [z/2*(r[i]-q_geometric[i])-(r[i]-1) for i in range(len(x))]

    gamma = ()

    for i in range(len(x)):

        # Calculating first term in residual contribution for the component
        interaction_term_1 = numpy.dot(theta_interaction, numpy.transpose(tau)[i])

        # Calculating second term in residual contribution for the component
        interaction_term_2 = 0

        for j in range(len(x)):

            interaction_term_2_sum = 0
            for k in range(len(x)):
                interaction_term_2_sum += theta_interaction[k]*tau[k][j]

            interaction_term_2 += theta_interaction[j]*tau[i][j]/interaction_term_2_sum

        gamma += ((numpy.exp(
            numpy.log(phi[i]/x[i])
            + z/2 * q_geometric[i] * numpy.log(theta_geometric[i] / phi[i])
            + l[i]
            - phi[i] / x[i] * numpy.sum(numpy.multiply(l, x))
            - q_interaction[i] * (
                + numpy.log(interaction_term_1)
                + interaction_term_2
                - 1)
         )
        ),)

    return gamma
