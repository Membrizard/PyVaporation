from pyvaporation.mixtures import Mixtures, get_partial_pressures, Composition
from pyvaporation.mixtures import VLEPoints, fit_vle

points = VLEPoints.from_csv(path="./tests/VLE_data/ternary/H2O_MeOH_EtOH.csv")

params = fit_vle(data=points)

print(params)
