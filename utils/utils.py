import typing

import attr
import numpy

from ..component import Component

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
    permeate_temperature: float
    feed_amount: float
    feed_composition: Composition
    isothermal: bool = True


@attr.s(auto_attribs=True)
class AntoineConstants:
    a: float = attr.ib(converter=lambda value: float(value))  # type: ignore
    b: float = attr.ib(converter=lambda value: float(value))  # type: ignore
    c: float = attr.ib(converter=lambda value: float(value))  # type: ignore


@attr.s(auto_attribs=True)
class NRTLParameters:
    g12: float
    g21: float
    alpha: float


@attr.s(auto_attribs=True)
class Mixture:
    components: typing.List[Component]
    nrtl_params: NRTLParameters

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

    # H2O/Ethanol
    @classmethod
    def h2o_etoh(cls) -> 'Mixture':
        return cls(
            components=[Component.h2o(), Component.etoh()],
            nrtl_params=NRTLParameters(
                g12=-633,
                g21=5823,
                alpha=0.3,
            )
        )

    # H2O/Isopropanol TODO: update parameters
    @classmethod
    def h2o_ipoh(cls) -> 'Mixture':
        return cls(
            components=[Component.h2o(), Component.ipoh()],
            nrtl_params=NRTLParameters(
                g12=0,
                g21=0,
                alpha=0.686,
            )
        )

    # EtOH/ETBE
    @classmethod
    def etoh_etbe(cls) -> 'Mixture':
        return cls(
            components=[Component.etoh(), Component.etbe()],
            nrtl_params=NRTLParameters(
                g12=726.34,
                g21=318.22,
                alpha=0.3,
            )
        )

    # MeOH/Toluene
    @classmethod
    def meoh_toluene(cls) -> 'Mixture':
        return cls(
            components=[Component.meoh(), Component.toluene()],
            nrtl_params=NRTLParameters(
                g12=3715.5266,
                g21=3085.3254,
                alpha=0.3,
            )
        )

    # MeOH/Methyl-tertButhyl ether
    @classmethod
    def meoh_mtbe(cls) -> 'Mixture':
        return cls(
            components=[Component.meoh(), Component.mtbe()],
            nrtl_params=NRTLParameters(
                g12=2025.3,
                g21=2133.5,
                alpha=0.6,
            )
        )

    # MeOH/DimethoxyEthane
    @classmethod
    def meoh_dme(cls) -> 'Mixture':
        return cls(
            components=[Component.meoh(), Component.dme()],
            nrtl_params=NRTLParameters(
                g12=782.0202,
                g21=-229.0431,
                alpha=0.2982,
            )
        )

    # MeOH/DimethylCarbonate
    @classmethod
    def meoh_dmc(cls) -> 'Mixture':
        return cls(
            components=[Component.meoh(), Component.dmc()],
            nrtl_params=NRTLParameters(
                g12=1429.444,
                g21=2641.108,
                alpha=0.2,
            )
        )

    # MeOH/Cyclohexane
    @classmethod
    def meoh_cyclohexane(cls) -> 'Mixture':
        return cls(
            components=[Component.meoh(), Component.cyclohexane()],
            nrtl_params=NRTLParameters(
                g12=6415.36,
                g21=5714,
                alpha=0.4199,
            )
        )


# Experiments for Ideal modelling
@attr.s(auto_attribs=True)
class IdealExperiment:
    temperature: float
    # Permeance in kg*mcm/(m2*h*kPa)
    permeance: float
    component: Component
    activation_energy: typing.Optional[float] = None


# Experiments for Non-ideal modelling
@attr.s(auto_attribs=True)
class NonIdealExperiment:
    temperature: float
    # Permeance in kg*mcm/(m2*h*kPa)
    overall_flux: typing.List[float]
    component_fluxes: typing.List[typing.List[float]]
    compositions: typing.List[float]
    permeate_temperature: typing.Optional[float] = 0
    mixture: Mixture


@attr.s(auto_attribs=True)
class Membrane:
    ideal_experiments: typing.List[IdealExperiment]
    non_ideal_experiments: typing.List[NonIdealExperiment]

    @property
    def ideal_experiments_names(self) -> typing.List[str]:
        return [ie.component.name for ie in self.ideal_experiments]

    # Get all the penetrants the membrane was tested for
    def get_known_penetrants(self) -> typing.List[Component]:
        return numpy.unique([
            ideal_experiment.component for ideal_experiment in self.ideal_experiments
        ])

    # Picking only Experiments related to the chosen component
    def get_penetrant_data(self, component) -> typing.List[IdealExperiment]:
        return list(
            filter(
                lambda value: value.component.name in self.ideal_experiments_names,
                self.ideal_experiments
            )
        )

    # Calculate an apparent activation energy of permeation
    def calculate_activation_energy(self, component, rewrite: bool = True, plot: bool = False) -> float:
        # Validation of the Ideal Experiments
        if len(self.get_penetrant_data()) >= 2:
            # Get all the temperature values corresponding to the permeances for a given penetrant convert to 1/T
            x = numpy.power([ideal_experiment.temperature for ideal_experiment in
                             self.get_penetrant_data(component)], -1)
            # Get all the permeance values for a given Penetrant convert to Ln(permeance)
            y = numpy.log([ideal_experiment.permeance for ideal_experiment in
                           self.get_penetrant_data(component)])
            # Converting Ln(permeance) to the matrix equation form
            A = numpy.vstack([x, numpy.ones(len(x))]).T
            # Calculation of Least-squares Linear Fit of Ln(Permeance) versus 1/T
            activation_energy, c = numpy.linalg.lstsq(A, y, rcond=-1)[0] * R
            if plot:
                # TODO: Plotting the graph Ln(Permeance) versus 1/T
                pass
            if rewrite:
                # Rewriting the corresponding activation energy values in Ideal Experiments of the Membrane,
                # including in Config
                return activation_energy
            else:
                return activation_energy
        else:
            print("At Least Two points are required for the calculation of Apparent Activation Energy")

    def get_permeance(self, temperature, component) -> float:
        return 0

    def get_ideal_selectivity(self, temperature, component1, component2) -> float:
        return 0


@attr.s(auto_attribs=True)
class Pervaporation:
    membrane: Membrane
    mixture: Mixture
    conditions: Conditions
    ideal: bool = True

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

    def model_ideal_diffusion_curve(self):
        pass

    def model_ideal_process(self):
        pass

    def model_non_ideal_process(self):
        pass
