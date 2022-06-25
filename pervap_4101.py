from pathlib import Path
from matplotlib import pyplot as plt

from conditions import Conditions
from mixtures import Mixtures, Composition, CompositionType, get_nrtl_partial_pressures
from pervaporation import Pervaporation
from membrane import Membrane
from config import Config
from permeance import Permeance


config = Config(
    source_path=Path("tests/data/Pervap_4101"),
)

membrane = Membrane.load(config)

pv = Pervaporation(
    membrane=membrane,
    mixture=Mixtures.H2O_EtOH,
)

con = Conditions(
    membrane_area=0.017,
    initial_feed_temperature=368.15,
    initial_feed_amount=1.5,
    initial_feed_composition=Composition(p=0.1, type=CompositionType.weight),
    permeate_pressure=0,
)

model = pv.non_ideal_isothermal_process(
    conditions=con,
    diffusion_curve_set=membrane.diffusion_curve_sets[0],
    initial_permeances=(
        Permeance(0.0143),
        Permeance(0.00000632),
    ),
    number_of_steps=50,
    delta_hours=0.2,
)

experiment_time_hours = [
    0,
    0.5163,
    1.0148,
    1.5133,
    2.0118,
    3.0087,
    4.0057,
    5.0173,
    6.0143,
    7.0112,
    8.0082,
    9.0052,
    10,
]

experiment_water_fraction = [
    0.1,
    0.0946,
    0.0938,
    0.0907,
    0.0887,
    0.0842,
    0.0791,
    0.0730,
    0.0702,
    0.0669,
    0.0625,
    0.0588,
    0.0573,
]

pressures = get_nrtl_partial_pressures(
    368.15, Mixtures.H2O_EtOH, Composition(p=0.1, type=CompositionType.weight)
)

experiment_water_flux = [
    0.4863,
    0.5209,
    0.5301,
    0.5865,
    0.4517,
    0.4736,
    0.4353,
    0.3789,
    0.3770,
    0.3370,
    0.2987,
    0.2823,
    0.2823,
]
print(experiment_water_flux[0] / pressures[0])


plt.plot(experiment_time_hours, experiment_water_fraction)
plt.plot(model.time, [c.first for c in model.feed_composition])
plt.legend(["experiment", "model"])
plt.show()

plt.plot(experiment_time_hours, experiment_water_flux)
plt.plot(model.time, [p[0] for p in model.partial_fluxes])
plt.legend(["experiment", "model"])
plt.show()