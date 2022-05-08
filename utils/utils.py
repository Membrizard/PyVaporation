import typing

import attr
import numpy

R = 8.314462


@attr.s(auto_attribs=True)
class Composition:
    # TODO: from 2 to n
    p: float = attr.ib(validator=lambda value: 0 <= value <= 1)  # type: ignore

    def __getitem__(self, item: int):
        if item == 0:
            return self.p
        elif item == 1:
            return 1 - self.p
        else:
            raise ValueError("Index %s out of range" % item)


@attr.s(auto_attribs=True)
class Conditions:
    membrane_area: float
    temperature: float
    permeate_pressure: float
    feed_amount: float
    feed_composition: Composition
    isothermal: bool = True


@attr.s(auto_attribs=True)
class AntoineConstants:
    a: float
    b: float
    c: float


@attr.s(auto_attribs=True)
class NRTLParameters:
    g12: float
    g21: float
    alpha: float


@attr.s(auto_attribs=True)
class Component:
    name: str
    antoine_constants: AntoineConstants

    def get_antoine_pressure(self, temperature: float) -> float:
        return 100 * 10 ** (
            self.antoine_constants.a
            - self.antoine_constants.b / (temperature + self.antoine_constants.c)
        )

    def get_vaporization_heat(self, temperature: float) -> float:
        return (
            (temperature / (temperature + self.antoine_constants.c)) ** 2
            * R
            * self.antoine_constants.b
            * numpy.log(10)
        )


@attr.s(auto_attribs=True)
class Mixture:
    components: typing.List[Component]
    nrtl_params: NRTLParameters

    def get_nrtl_partial_pressures(self, temperature: float, composition: Composition):

        g_scaled = numpy.array(
            [
                self.nrtl_params.g12 / (R * temperature),
                self.nrtl_params.g21 / (R * temperature),
            ]
        )

        g_exp = numpy.exp(-g_scaled * self.nrtl_params.alpha)

        activity_coefficients = [
            numpy.exp(
                (composition[1] ** 2)
                * (
                    g_scaled[1]
                    * (g_exp[1] / (composition[0] + composition[1] * g_exp[1])) ** 2
                    + g_scaled[0]
                    * g_exp[0]
                    / (composition[1] + composition[0] * g_exp[0]) ** 2
                )
            ),
            numpy.exp(
                (composition[0] ** 2)
                * (
                    g_scaled[0]
                    * (g_exp[0] / (composition[1] + composition[0] * g_exp[0])) ** 2
                    + g_scaled[1]
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
class Experiment:
    temperature: float
    permeance: float
    component: Component
    activation_energy: typing.Optional[float] = None


@attr.s(auto_attribs=True)
class Membrane:
    experiments: typing.List[Experiment]


    def calculate_activation_energy(self) -> float:
        return 0


@attr.s(auto_attribs=True)
class Pervaporation:
    membrane: Membrane
    mixture: Mixture
    conditions: Conditions

    def get_component_flux(self):
        pass

    def get_overall_flux(self):
        pass

    def get_separation_factor(self):
        pass

    def get_psi(self):
        pass

    def get_real_selectivity(self):
        pass

    def model_process(self):
        pass
