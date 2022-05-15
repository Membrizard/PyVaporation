import attr
import numpy

from ..membrane import Membrane
from ..mixture import Mixture
from ..conditions import Conditions


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
