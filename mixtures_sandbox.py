from pyvaporation.mixtures import Mixtures, get_partial_pressures, Composition
from pyvaporation.components import Components

mixture = Mixtures.H2O_EtOH



ps = get_partial_pressures(temperature=293,
                           mixture=mixture,
                           composition=Composition(p=0.5, type="molar"),
                           )


print(ps)
