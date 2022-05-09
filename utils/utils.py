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
    permeate_temperature: float
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
    molecular_weight: float
    antoine_constants: AntoineConstants
    # Calculation of saturated pressure in kPa at a given temperature in K using Antoine equation (by the basis of 10)
    def get_antoine_pressure(self, temperature: float) -> float:
        return 10 ** (
            self.antoine_constants.a
            - self.antoine_constants.b / (temperature + self.antoine_constants.c)
        )
    # Calculation of Vaporisation heat in kJ/mol using Clapeyron-Clausius equation
    def get_vaporisation_heat(self, temperature: float) -> float:
        return (
            (temperature / (temperature + self.antoine_constants.c)) ** 2
            * R
            * self.antoine_constants.b
            * numpy.log(10)
        )
    
    # Water
    @classmethod
    def h2o(cls) -> 'Component':
        return cls(
            name='H2O',
            molecular_weight=18.02,
            antoine_constants=AntoineConstants(
                a=7.20389,
                b=1733.96,
                c=-39.485,
            )
        )
   # Methanol
    @classmethod
    def meoh(cls) -> 'Component':
        return cls(
            name='MeOH',
            molecular_weight=32.04,
            antoine_constants=AntoineConstants(
                a=9.22010,
                b=-1590.15535,
                c=-32.77001,
            )
        )
    # Ethanol
    @classmethod
    def etoh(cls) -> 'Component':
        return cls(
            name='EtOH',
            molecular_weight=46.07,
            antoine_constants=AntoineConstants(
                a=7.24677,
                b=1598.673,
                c=-46.424,
            )
        )
    # Isopropanol
    @classmethod
    def ipoh(cls) -> 'Component':
        return cls(
            name='IPOH',
            molecular_weight=60.10,
            antoine_constants=AntoineConstants(
                a=6.861,
                b=1357.427,
                c=-75.814,
            )
        )
   # Methyl-tertbuthyl ether
    @classmethod
    def mtbe(cls) -> 'Component':
        return cls(
            name='MTBE',
            molecular_weight=88.15,
            antoine_constants=AntoineConstants(
                a=8.05093,
                b=-1139.816725,
                c=-46.424,
            )
        )
    # Ethyl-tertbuthyl ether
    @classmethod
    def etbe(cls) -> 'Component':
        return cls(
            name='ETBE',
            molecular_weight=102.17,
            antoine_constants=AntoineConstants(
                a=6.0703,
                b=-1206.874,
                c=-49.19,
            )
        )
   # DimethoxyEthane
    @classmethod
    def dme(cls) -> 'Component':
        return cls(
            name='DME',
            molecular_weight=90.12,
            antoine_constants=AntoineConstants(
                a=5.83775,
                b=-1260.52,
                c=-37.322,
            )
        )
    # DimethylCarbonate
    @classmethod
    def dmc(cls) -> 'Component':
        return cls(
            name='DMC',
            molecular_weight=90.08,
            antoine_constants=AntoineConstants(
                a=5.78894,
                b=-1049.375,
                c=-85.977,
            )
        )
    # Cyclohexane
    @classmethod
    def cyclohexane(cls) -> 'Component':
        return cls(
            name='Cyclohexane',
            molecular_weight=84.16,
            antoine_constants=AntoineConstants(
                a=5.97636,
                b=1206.47,
                c=-49.864,
            )
        )
    # Benzene
    @classmethod
    def benzene(cls) -> 'Component':
        return cls(
            name='Benzene',
            molecular_weight=78.11,
            antoine_constants=AntoineConstants(
                a=6.00477,
                b=1196.76,
                c=-53.839,
            )
        )
    # Toluene
    @classmethod
    def toluene(cls) -> 'Component':
        return cls(
            name='Toluene',
            molecular_weight=92.14,
            antoine_constants=AntoineConstants(
                a=6.0854,
                b=1348.77,
                c=-53.024,
            )
        )
    # p-Xylene
    @classmethod
    def p_xylene(cls) -> 'Component':
        return cls(
            name='p-Xylene',
            molecular_weight=106.16,
            antoine_constants=AntoineConstants(
                a=8.11543,
                b=1453.43,
                c=-57.7,
            )
        )

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
    def h2o_etoh(cls)-> 'Mixture':
        return cls(
            components = [Component.h2o(),Component.etoh()],
            nrtl_params = NRTLParameters (
                g12 = -633,
                g21 = 5823,
                alpha = 0.3,
            )
        )
    # H2O/Isopropanol TODO: update parameters
    @classmethod
    def h2o_ipoh(cls)-> 'Mixture':
        return cls(
            components = [Component.h2o(),Component.ipoh()],
            nrtl_params = NRTLParameters (
                g12 = 0,
                g21 = 0,
                alpha = 0.686,
            )
        )
    # EtOH/ETBE
    @classmethod
    def etoh_etbe(cls)-> 'Mixture':
        return cls(
            components = [Component.etoh(),Component.etbe()],
            nrtl_params = NRTLParameters (
                g12 = 726.34,
                g21 = 318.22,
                alpha = 0.3,
            )
        )
    # MeOH/Toluene
    @classmethod
    def meoh_toluene(cls)-> 'Mixture':
        return cls(
            components = [Component.meoh(),Component.toluene()],
            nrtl_params = NRTLParameters (
                g12 = 3715.5266,
                g21 = 3085.3254,
                alpha = 0.3,
            )
        )
    # MeOH/Methyl-tertButhyl ether
    @classmethod
    def meoh_mtbe(cls)-> 'Mixture':
        return cls(
            components = [Component.meoh(),Component.mtbe()],
            nrtl_params = NRTLParameters (
                g12 = 2025.3,
                g21 = 2133.5,
                alpha = 0.6,
            )
        )
    # MeOH/DimethoxyEthane
    @classmethod
    def meoh_dme(cls)-> 'Mixture':
        return cls(
            components = [Component.meoh(),Component.dme()],
            nrtl_params = NRTLParameters (
                g12 = 782.0202,
                g21 = -229.0431,
                alpha = 0.2982,
            )
        )
    # MeOH/DimethylCarbonate
    @classmethod
    def meoh_dmc(cls)-> 'Mixture':
        return cls(
            components = [Component.meoh(),Component.dmc()],
            nrtl_params = NRTLParameters (
                g12 = 1429.444,
                g21 = 2641.108,
                alpha = 0.2,
            )
        )
    # MeOH/Cyclohexane
    @classmethod
    def meoh_cyclohexane(cls)-> 'Mixture':
        return cls(
            components = [Component.meoh(),Component.cyclohexane()],
            nrtl_params = NRTLParameters (
                g12 = 6415.36,
                g21 = 5714,
                alpha = 0.4199,
            )
        )
# Experiments for Ideal modelling
@attr.s(auto_attribs=True)
class Ideal_Experiment:
    temperature: float
    # Permeance in kg*mcm/(m2*h*kPa)
    permeance: float
    component: Component
    activation_energy: typing.Optional[float] = None

# Experiments for Non-ideal modelling
@attr.s(auto_attribs=True)
class NonIdeal_Experiment:
    temperature: float
    # Permeance in kg*mcm/(m2*h*kPa)
    overall_flux: typing.List[float]
    component_fluxes:typing.List[typing.List[float]]
    compositions: typing.List[float]
    mixture: Mixture

@attr.s(auto_attribs=True)
class Membrane:
    ideal_experiments: typing.List[Ideal_Experiment]
    nonideal_experiments: typing.List[NonIdeal_Experiment]

# Get all the penetrants the membrane was tested for
    def get_known_penetrants(self)-> typing.List[Component]:
        return numpy.unique([
            Ideal_Experiment.component() for Ideal_Experiment in self.ideal_experiments()
        ])
 # Picking only Experiments related to the chosen component
    def get_penetrant_data(self, component)-> typing.List[Ideal_Experiment]:
        def check(component)-> bool:
           if component in [Ideal_Experiment.component() for Ideal_Experiment in self.ideal_experiments()]:
            return True
           else:
            return False
        return filter(check(component),self.ideal_experiments)

   # Calculate an apparent activation energy of permeation     
    def calculate_activation_energy(self,component) -> float:
        
        self.get_penetrant_data(component)
        # Calculation of Least-squares Linear Fit of Ln(Permeance) versus 1/T
        return 0


    def get_permeance(self,temperature,component)-> float:
        return 0
    def get_ideal_selectivity(self,temperature,component1,component2)-> float:
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

    def model_nonideal_process(self):
        pass
