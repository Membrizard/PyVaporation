from mixtures import Mixtures, Composition, CompositionType, get_nrtl_partial_pressures
from components import Components
from diffusion_curve import DiffusionCurve, DiffusionCurveSet
from membrane import Membrane
from mixtures import Composition, CompositionType, Mixtures
from permeance import Permeance
from optimizer import Measurements, find_best_fit
from pervaporation import Pervaporation
from conditions import Conditions
from config import Config
from pathlib import Path
import matplotlib.pyplot as plt
import numpy

config = Config(
        source_path=Path("tests/data/Pervap_2510"),
    )

membrane = Membrane.load(config)

pv = Pervaporation(
        membrane=membrane,
        mixture=Mixtures.H2O_iPOH,
    )

con = Conditions(
        membrane_area=0.017,
        initial_feed_temperature=368.15,
        initial_feed_amount=1.5,
        initial_feed_composition=Composition(p=0.15, type=CompositionType.weight),
        permeate_pressure=0,
    )

# model = pv.non_ideal_isothermal_process(
#         conditions=con,
#         diffusion_curve_set=membrane.diffusion_curve_sets[0],
#         number_of_steps=50,
#         delta_hours=0.2,
#     )

measurements = Measurements.from_diffusion_curves_first(membrane.diffusion_curve_sets[0])
fit = find_best_fit(measurements)
fit.plot()
fit.plot(concentration=(0.1, 0.2))
fit.plot(temperature=368.15)
fit.plot(concentration=(0.1, 0.2), temperature=368.15)
fit.plot(measurements)
fit.plot(measurements, temperature=363.15)
fit.plot(measurements, concentration=(0, 1))
fit.plot(measurements, temperature=368.15)


