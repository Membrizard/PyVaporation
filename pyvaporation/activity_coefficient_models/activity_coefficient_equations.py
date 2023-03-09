import numpy
import typing
from ..utils import NRTLParameters, R, UNIQUACConstants


class ActivityCoefficientModel:
    """
    Class to represent types of activity coefficient models used for calculation of partial pressures
    """

    NRTL: str = "NRTL"
    UNIQUAC: str = "UNIQUAC"


def uniquac_activity_coefficient_equation(
    # uniquac_constants: typing.List[UNIQUACConstants],
    # binary_interaction_parameters_matrix: typing.List[typing.List[float]],
    r: typing.List[float],
    q_geometric: typing.List[float],
    q_interaction: typing.List[float],
    x: typing.List[float],
    tau: typing.List[typing.List[float]],
    z: int,
) -> typing.Tuple:
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
    if len(x) < 2:
        raise ValueError(
            "There should be at least two compositions two calculate\
         activity coefficients using UNIQUAC model"
        )

    phi_sum = numpy.dot(r, x)
    phi = [r[i] * x[i] / phi_sum for i in range(len(x))]

    theta_geometric_sum = numpy.dot(q_geometric, x)
    theta_geometric = [
        q_geometric[i] * x[i] / theta_geometric_sum for i in range(len(x))
    ]

    theta_interaction_sum = numpy.dot(q_interaction, x)
    theta_interaction = [
        q_interaction[i] * x[i] / theta_interaction_sum for i in range(len(x))
    ]

    l = [z / 2 * (r[i] - q_geometric[i]) - (r[i] - 1) for i in range(len(x))]

    gamma = ()

    for i in range(len(x)):

        # Calculating first term in residual contribution for the component
        interaction_term_1 = numpy.dot(theta_interaction, numpy.transpose(tau)[i])

        # Calculating second term in residual contribution for the component
        interaction_term_2 = 0

        for j in range(len(x)):

            interaction_term_2_sum = 0
            for k in range(len(x)):
                interaction_term_2_sum += theta_interaction[k] * tau[k][j]

            interaction_term_2 += (
                theta_interaction[j] * tau[i][j] / interaction_term_2_sum
            )

        gamma += (
            (
                numpy.exp(
                    numpy.log(phi[i] / x[i])
                    + z / 2 * q_geometric[i] * numpy.log(theta_geometric[i] / phi[i])
                    + l[i]
                    - phi[i] / x[i] * numpy.sum(numpy.multiply(l, x))
                    - q_interaction[i]
                    * (+numpy.log(interaction_term_1) + interaction_term_2 - 1)
                )
            ),
        )

    return gamma


def nrtl_activity_coefficient_equation(
    nrtl_params: NRTLParameters,
    temperature: float,
    x: typing.List[float],
) -> typing.Tuple:
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
    if len(x) != 2:
        raise ValueError("NRTL model supports calculations only for binary mixtures.")

    tau = numpy.array(
        [
            (nrtl_params.a12 + nrtl_params.g12 / (R * temperature)),
            (nrtl_params.a21 + nrtl_params.g21 / (R * temperature)),
        ]
    )
    if nrtl_params.alpha21 is None:
        alphas = nrtl_params.alpha12
    else:
        alphas = [nrtl_params.alpha12, nrtl_params.alpha21]

    g_exp = numpy.exp(numpy.multiply(-tau, alphas))

    gamma = (
        numpy.exp(
            (x[1] ** 2)
            * (
                tau[1] * (g_exp[1] / (x[0] + x[1] * g_exp[1])) ** 2
                + tau[0] * g_exp[0] / (x[1] + x[0] * g_exp[0]) ** 2
            )
        ),
        numpy.exp(
            (x[0] ** 2)
            * (
                tau[0] * (g_exp[0] / (x[1] + x[0] * g_exp[0])) ** 2
                + tau[1] * g_exp[1] / (x[0] + x[1] * g_exp[1]) ** 2
            )
        ),
    )

    return gamma
