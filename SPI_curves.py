from membrane import Membrane
from mixtures import Mixtures, Composition, CompositionType
from permeance import Permeance, Units
from pervaporation import Pervaporation
from diffusion_curve import DiffusionCurveSet, DiffusionCurve
from conditions import Conditions
import matplotlib.pyplot as plt
import numpy
from experiments import IdealExperiment, IdealExperiments
from components import Components
from optimizer import Measurements, find_best_fit

compositions = [0.040, 0.037, 0.030, 0.026, 0.024, 0.021]

flux_h2o = [
    0.062062,
    0.062062,
    0.03988709016,
    0.03720394737,
    0.004707604167,
    0.005613600498,
]

flux_etoh = [
    0.006138,
    0.006138,
    0.003944877049,
    0.00413377193,
    0.0007663541667,
    0.0009520141196,
]

experimental_time = [0, 1.03, 2.02, 2.97, 4.00, 5.02]

curve = DiffusionCurve(
    mixture=Mixtures.H2O_EtOH,
    membrane_name="Romakon-Al2",
    feed_temperature=319.65,
    feed_compositions=[Composition(c, CompositionType.weight) for c in compositions],
    partial_fluxes=[(flux_h2o[i], flux_etoh[i]) for i in range(len(flux_h2o))],
)

curve_set = DiffusionCurveSet(name_of_the_set=" ", diffusion_curves=[curve])

romakon_al2 = Membrane(name="Romakon Al2", diffusion_curve_sets=[curve_set])

pervaporation = Pervaporation(romakon_al2, Mixtures.H2O_EtOH)

conditions = Conditions(
    membrane_area=0.0048,
    initial_feed_composition=Composition(p=0.04, type=CompositionType.weight),
    initial_feed_amount=0.047,
    initial_feed_temperature=319.65,
)

modelled_experiment = pervaporation.non_ideal_isothermal_process(
    conditions=conditions,
    diffusion_curve_set=curve_set,
    number_of_steps=50,
    delta_hours=0.1,
)
c = [c.first for c in modelled_experiment.feed_composition]
c.pop(-1)
plt.plot(experimental_time, compositions)
plt.plot(modelled_experiment.time, c)
plt.show()
