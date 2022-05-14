import attr

from ..membrane import Membrane
from ..mixture import Mixture
from ..conditions import Conditions


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
