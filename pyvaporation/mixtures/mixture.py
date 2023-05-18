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
            nrtl_params=self.nrtl_params,
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
    nrtl_params: typing.Optional[NRTLParameters] = None
    uniquac_params: typing.Optional[UNIQUACParameters] = None

    def __len__(self):
        return len(self.components)

    @property
    def first_component(self):
        return self.components[0]

    @property
    def second_component(self):
        return self.components[1]

    def to_binary_mixture(self) -> "BinaryMixture":
        if len(self.components) != 2:
            raise ValueError(
                "Mixture should have 2 Components to be converted to a BinaryMixture."
            )
        return BinaryMixture(
            name=self.name,
            first_component=self.components[0],
            second_component=self.components[1],
            nrtl_params=self.nrtl_params,
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
                raise ValueError(f"{name} is not in Mixture Components: {names}")


class CompositionType:
    """
    A class to describe type of the composition
    """

    molar: str = "molar"
    weight: str = "weight"


@attr.s(auto_attribs=True)
class Composition:
    """
    A class to represent composition of the mixtures
    :param p
    """

    p: typing.Union[typing.List, float]
    type: str

    def __attrs_post_init__(self):
        p = self.p

        if isinstance(p, float) or isinstance(p, int):
            self.p = [p, 1 - p]

        if isinstance(p, list):
            for value in p:
                if not isinstance(value, float):
                    raise ValueError(f"{value} is not float.")

        for value in self.p:
            if not 0 <= value <= 1:
                raise ValueError(f"Given {value} value is not in [0, 1] range")
        if abs(1 - sum(self.p)) > 1e-10:
            raise ValueError("The Sum of the component fractions should be equal to 1")

    def __len__(self):
        return len(self.p)

    def __getitem__(self, item):
        return self.p[item]

    def __setitem__(self, key, value):
        self.p[key] = value

    @property
    def first(self) -> float:
        """
        Returns fraction of the first test_components
        """
        return self.p[0]

    @property
    def second(self) -> float:
        """
        Returns fraction of the second test_components
        """
        return self.p[1]

    def to_molar(self, mixture: typing.Union[Mixture, BinaryMixture]) -> "Composition":
        """
        Converts Composition to molar %
        """
        if isinstance(mixture, BinaryMixture):
            mixture = mixture.to_mixture()

        if self.type == CompositionType.molar:
            return self
        else:
            inv_molecular_weights = [
                1 / component.molecular_weight for component in mixture.components
            ]
            p = numpy.multiply(self.p, inv_molecular_weights)
            sum_p = sum(p)
            p = [value / sum_p for value in p]

        return Composition(p=p, type=CompositionType.molar)

    def to_weight(self, mixture: typing.Union[Mixture, BinaryMixture]) -> "Composition":
        """
        Converts Composition to weight %
        """
        if isinstance(mixture, BinaryMixture):
            mixture = mixture.to_mixture()

        if self.type == CompositionType.weight:
            return self
        else:
            molecular_weights = [
                component.molecular_weight for component in mixture.components
            ]
            p = numpy.multiply(self.p, molecular_weights)
            sum_p = sum(p)
            p = [value / sum_p for value in p]

            return Composition(p=p, type=CompositionType.weight)


def get_partial_pressures(
    temperature: float,
    mixture: typing.Union[BinaryMixture, Mixture],
    composition: Composition,
    calculation_type: str = ActivityCoefficientModel.NRTL,
) -> typing.Tuple[float]:
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
    pc_pressures = [
        component.get_vapor_pressure(temperature) for component in mixture.components
    ]
    r_pressures = numpy.multiply(pc_pressures, activity_coefficients)
    r_pressures = tuple(numpy.multiply(r_pressures, composition.p))
    return r_pressures


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

        for i in range(len(composition)):
            if composition[i] == 0:
                composition[i] = 1e-10

        if mixture.uniquac_params is None:
            raise ValueError(
                "UNIQUAC Parameters must be specified for this type of calculation"
            )
        component_consts = []
        for component in mixture.components:
            if component.uniquac_constants is None:
                raise ValueError(
                    "UNIQUAC Constants for all Components must be specified for this type of calculation"
                )
            component_consts.append(component.uniquac_constants)

        activity_coefficients = uniquac_activity_coefficient_equation(
            r=[component.r for component in component_consts],
            q_geometric=[component.q_geometric for component in component_consts],
            q_interaction=[component.q_interaction for component in component_consts],
            x=composition.p,
            tau=mixture.uniquac_tau_matrix(temperature=temperature),
            z=mixture.uniquac_params.z,
        )

        return activity_coefficients
