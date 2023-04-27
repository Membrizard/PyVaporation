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


@attr.s(auto_attribs=True, frozen=True)
class BinaryMixture:
    """
    A class to represent binary mixtures
    TODO lookup attr frozen how to apply
    """

    name: str
    first_component: Component
    second_component: Component
    nrtl_params: typing.Optional[NRTLParameters] = None
    uniquac_params: typing.Optional[UNIQUACParameters] = None

    def to_mixture(self) -> "Mixture":
        return Mixture(
            name=self.name,
            components=[self.first_component, self.second_component],
            # nrtl_params=self.nrtl_params,
            uniquac_params=self.uniquac_params,
        )

    def __len__(self):
        return 2

    def uniquac_tau_matrix(self, temperature: float) -> numpy.ndarray:
        interaction_matrix = numpy.ones(shape=(2, 2))

        parameters = self.uniquac_params.binary_parameters[0]

        ij_parameter = parameters.ij_parameter
        ji_parameter = parameters.ji_parameter
        ij_tau = numpy.exp(
            -(ij_parameter[0] + ij_parameter[1] / temperature) / temperature
        )
        ji_tau = numpy.exp(
            -(ji_parameter[0] + ji_parameter[1] / temperature) / temperature
        )
        interaction_matrix[0][1] = ij_tau
        interaction_matrix[1][0] = ji_tau

        return interaction_matrix

    def __attrs_post_init__(self):
        if self.nrtl_params is None and self.uniquac_params is None:
            raise ValueError(
                "Component Interaction parameters are required to create a mixture!"
            )

    def change_components(
        self,
        first_component: typing.Optional[Component] = None,
        second_component: typing.Optional[Component] = None,
        uniquac_params: typing.Optional[UNIQUACParameters] = None,
        nrtl_params: typing.Optional[NRTLParameters] = None,
    ) -> "BinaryMixture":

        if first_component == second_component:
            raise ValueError("BinaryMixture should contain 2 different components")

        components = [self.first_component, self.second_component]
        if first_component is None and second_component is None:
            return self

        if first_component in components or second_component in components:
            if first_component != second_component:
                return self

        if first_component is not None or second_component is not None:
            if uniquac_params is None and nrtl_params is None:
                raise ValueError(
                    "Interaction parameters must be updated when components are changed."
                )
            return BinaryMixture(
                name=self.name,
                first_component=first_component,
                second_component=second_component,
                nrtl_params=nrtl_params,
                uniquac_params=uniquac_params,
            )


@attr.s(auto_attribs=True, frozen=True)
class Mixture:
    """
    A class to represent multicomponent mixtures
    """

    name: str
    components: typing.List[Component]
    # nrtl_params: typing.Optional[NRTLParameters] = None
    uniquac_params: UNIQUACParameters

    def __len__(self):
        return len(self.components)

    def to_binary_mixture(self) -> "BinaryMixture":
        if len(self.components) != 2:
            raise ValueError(
                "Mixture should have 2 Components to be converted to a BinaryMixture."
            )
        return BinaryMixture(
            name=self.name,
            first_component=self.components[0],
            second_component=self.components[1],
            nrtl_params=None,
            uniquac_params=self.uniquac_params,
        )

    def uniquac_tau_matrix(self, temperature: float) -> numpy.ndarray:
        """
        The method calculates binary interaction parameters matrix at a given temperature (tau matrix)
        :param temperature: temperature in K
        :return: binary interaction parameters matrix
        """
        L = len(self)
        interaction_matrix = numpy.ones(shape=(L, L))
        parameters = self.uniquac_params.binary_parameters
        names = [component.name for component in self.components]
        parameter_coordinates = [
            (
                names.index(parameter.i_component_name),
                names.index(parameter.j_component_name),
            )
            for parameter in parameters
        ]

        for i in range(len(parameter_coordinates)):
            ij_parameter = parameters[i].ij_parameter
            ji_parameter = parameters[i].ji_parameter
            ij_tau = numpy.exp(
                -(ij_parameter[0] + ij_parameter[1] / temperature) / temperature
            )
            ji_tau = numpy.exp(
                -(ji_parameter[0] + ji_parameter[1] / temperature) / temperature
            )
            interaction_matrix[parameter_coordinates[i][0]][
                parameter_coordinates[i][1]
            ] = ij_tau
            interaction_matrix[parameter_coordinates[i][1]][
                parameter_coordinates[i][0]
            ] = ji_tau

        return interaction_matrix

    def __attrs_post_init__(self):

        if len(self) > 2 and (len(self) != len(self.uniquac_params)):
            raise ValueError(
                "Binary Interaction parameters matrix size does not correspond with the Mixture!"
            )

        names = [component.name for component in self.components]
        parameter_names = []

        for parameter in self.uniquac_params.binary_parameters:
            parameter_names.extend(
                [parameter.i_component_name, parameter.j_component_name]
            )

        for name in parameter_names:
            if not (name in names):
                raise ValueError(
                    f"{name} is not in Mixture Components: {names}"
                )


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

    def to_molar(self, mixture: BinaryMixture) -> "Composition":
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

    def to_weight(self, mixture: BinaryMixture) -> "Composition":
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
    mixture: typing.Union[BinaryMixture, Mixture],
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
    mixture: typing.Union[BinaryMixture, Mixture],
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

        if len(mixture) != 2:
            raise ValueError("NRTL model can be used only for Binary Mixtures")

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

        # binary_interaction_params = mixture.uniquac_params.binary_parameters_matrix
        # # tau[i][i] = 1 !!!!
        tau = mixture.uniquac_tau_matrix(temperature=temperature)

        # for i in range(len(binary_interaction_params)):
        #     for j in range(len(binary_interaction_params)):
        #         if binary_interaction_params[i][j] != 0:
        #             tau[i][j] = numpy.exp(
        #                 -(
        #                     binary_interaction_params[i][j][0]
        #                     + binary_interaction_params[i][j][1] / temperature
        #                 )
        #                 / temperature
        #             )

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
