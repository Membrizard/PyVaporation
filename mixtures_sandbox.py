from pyvaporation.mixtures import Mixtures, get_partial_pressures, Composition
from pyvaporation.components import Components

mixture = Mixtures.H2O_EtOH.to_binary_mixture()

print(mixture.uniquac_tau_matrix(293))

ps = get_partial_pressures(temperature=293,
                           mixture=mixture,
                           composition=Composition(p=0.5, type="molar"),
                           calculation_type="UNIQUAC")


print(ps)
