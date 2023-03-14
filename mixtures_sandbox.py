from pyvaporation.mixtures import Mixtures
from pyvaporation.components import Components

mixture = Mixtures.H2O_EtOH

print(mixture.uniquac_tau_matrix(293))
