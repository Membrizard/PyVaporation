import attr
import numpy

from enum import Enum

from ..mixture import Mixture

R = 8.314462


class CompositionType(Enum):
    molar: str = "molar"
    weight: str = "weight"


@attr.s(auto_attribs=True)
class Composition:
    p: float = attr.ib(validator=lambda value: 0 <= value <= 1)  # type: ignore
    type: CompositionType

    @property
    def first(self) -> float:
        return self.p

    @property
    def second(self) -> float:
        return 1 - self.p

    def to_molar(self, mixture: Mixture) -> "Composition":
        if self.type == CompositionType.molar:
            return self
        else:
            p = (self.p / mixture.first_component.molecular_weight) / (
                    self.p / mixture.first_component.molecular_weight
                    + (1 - self.p) / mixture.second_component.molecular_weight
            )
            return Composition(p=p, type=CompositionType('weight'))

    def to_weight(self, mixture: Mixture) -> "Composition":
        if self.type == CompositionType.weight:
            return self
        else:
            p = (mixture.first_component.molecular_weight * self.p) / (
                mixture.first_component.molecular_weight * self.p
                + mixture.second_component.molecular_weight * (1 - self.p)
            )
            return Composition(p=p, type=CompositionType('weight'))


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
class Membrane:
    ideal_experiments: typing.List[IdealExperiment]
    DiffusionCurve: typing.List[DiffusionCurve]

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
                # Plotting the graph Ln(Permeance) versus 1/T
                import matplotlib.pyplot as plt
                _ = plt.plot(x, y, 'o', label=component.name + 'Experimental Permeances', markersize=5)
                plt.xlabel("1/T, 1/K")
                plt.ylabel("Ln(Permeance)")
                _ = plt.plot(x, activation_energy * x + c, 'b', label='Fitted line')
                _ = plt.legend()
                plt.show()
                pass
            if rewrite:
                # Rewriting the corresponding activation energy values in Ideal Experiments of the Membrane
                for ideal_experiment in self.get_penetrant_data(component):
                    ideal_experiment.activation_energy = activation_energy
                # TODO: Rewrite in Config Membranes.yml or whatever

                return activation_energy
            else:
                return activation_energy
        else:
            print("At Least Two points are required for the calculation of Apparent Activation Energy")

    def get_permeance(self, temperature, component) -> float:
        # Definition of the corresponding temperatures list
        temperature_list = [ideal_experiment.temperature
                            for ideal_experiment in
                            self.get_penetrant_data(component)]
        # finding the index of the experiment with the closest available temperature
        index = min(range(len(temperature_list)), key=lambda i: abs(temperature_list[i] - temperature))
        # finding the Permeance in the experiment with the found index
        given_permeance = self.get_penetrant_data(component)[index].permeance
        # Trying to calculate the permeance at given temperature;
        # If activation is not specified it is being calculated using calculate_activation_energy function
        try:
            activation_energy = self.get_penetrant_data(component)[index].activation_energy
            return given_permeance * numpy.exp(-activation_energy / R * (1 / temperature - 1 / temperature_list[index]))
        except:
            # print('Provide Apparent Energy of Permeation for the component')
            activation_energy = self.calculate_activation_energy(component, rewrite=False, plot=False)
            return given_permeance * numpy.exp(-activation_energy / R * (1 / temperature - 1 / temperature_list[index]))

    # calculating the selectivity
    def get_ideal_selectivity(self, temperature, component1, component2) -> float:
        return self.get_permeance(temperature, component1) / self.get_permeance(temperature, component2)

# Diffusion curves for non-ideal modelling and output of the results
@attr.s(auto_attribs=True)
class DiffusionCurve:
    feed_temperature: float
    permeate_temperature: typing.Optional[float] = 0
    mixture: Mixture
    composition_range: typing.List[Composition]
    partial_fluxes: typing.List[typing.List[float]]
    total_flux: typing.List[float]
    separation_factor: typing.List[float]
    PSI: typing.List[float]

@attr.s(auto_attribs=True)
class PVProcess:


@attr.s(auto_attribs=True)
class Pervaporation:
    membrane: Membrane
    mixture: Mixture
    conditions: Conditions
    ideal: bool = True

    # Alexey please double check this
    def calculate_partial_fluxes(self, feed_temperature, permeate_temperature, composition, precision: float) -> float:
        # Calculating components' permeances at a given feed temperature:
        permeances = [self.membrane.get_permeance(feed_temperature, component) for component in self.mixture.components]
        # Defining function for saturation pressure calculation
        p_sat = lambda t, x: self.mixture.get_nrtl_partial_pressures(t, x)
        # Defining function for partial fluxes calculation from permeate composition
        partial_fluxes = lambda perm_comp: \
            numpy.matmul(permeances, numpy.substract(p_sat(feed_temperature, composition.to_mol()),
                                                     p_sat(permeate_temperature, perm_comp.to_mol)))
        # Defining function for permeate composition calculation
        permeate_composition = lambda fluxes: Composition(fluxes[0] / numpy.sum(fluxes))
        inital_fluxes = numpy.matmul(permeances, self.mixture.get_nrtl_partial_pressures(feed_temperature, composition))
        p_c = permeate_composition(inital_fluxes)
        d = 1
        # Precision of the permeate composition value - until the precision criteria is reached
        # That is definetly could be implemented in a better algorithmic way
        while d >= precision:
            p_c2 = permeate_composition(partial_fluxes(p_c))
            d = max(numpy.absolute(numpy.subtract(p_c2, p_c)))
            p_c = p_c2
        return partial_fluxes(p_c)

    # Calculate Permeate composition for at the given conditions
    def calculate_permeate_composition(self, feed_temperature, permeate_temperature, composition,
                                       precision: float) -> Composition:
        x = self.calculate_partial_fluxes(feed_temperature, permeate_temperature, composition, precision)
        return Composition(x[0] / numpy.sum(x))
    def calculate_separation_factor(self, feed_temperature, permeate_temperature, composition, precision: float):
        perm_comp = self.calculate_permeate_composition(feed_temperature,permeate_temperature,composition,precision)
        return (composition[1]/(1-composition[1]))/(perm_comp[1]/(1-perm_comp[1]))

    # Calculate Partial, Overall fluxes and other parameters as a function of composition in the given composition range
    def ideal_diffusion_curve(self, feed_temperature, permeate_temperature, composition_range, precision,
                              plot: bool = True) -> DiffusionCurve:
        diff_curv = DiffusionCurve()
        diff_curv.feed_temperature = feed_temperature
        diff_curv.mixture = self.mixture
        diff_curv.permeate_temperature = permeate_temperature
        diff_curv.composition_range_range = composition_range
        diff_curv.partial_fluxes = [self.calculate_partial_fluxes( feed_temperature, permeate_temperature, composition, precision) for
        composition in composition_range]
        diff_curv.separation_factor = [self.calculate_separation_factor( feed_temperature, permeate_temperature, composition, precision) for
        composition in composition_range]
        diff_curv.total_flux=[numpy.sum(diff_curv.partial_fluxes[i]) for i in range(len(diff_curv.partial_fluxes))]
        diff_curv.PSI = numpy.multiply(numpy.subtract(diff_curv.separation_factor, numpy.ones(0, len(composition_range))), diff_curv.total_flux)
        if plot:
            #TODO Plot the curve
            pass
        return diff_curv

    def model_ideal_process(self, conditions):
        pass

    def model_non_ideal_process(self, conditions):
        pass
