from optimizer import Measurements, fit, find_best_curve
import matplotlib.pyplot as plt
from components import Components
from conditions import Conditions
from experiments import IdealExperiment, IdealExperiments
from membrane import Membrane
from mixtures import Mixtures, Composition, CompositionType
from permeance import Permeance
from pervaporation import Pervaporation
from diffusion_curve import DiffusionCurveSet


experiment_h2o_1 = IdealExperiment(
    name="Romakon-PM102",
    temperature=323.15,
    component=Components.H2O,
    permeance=Permeance(0.036091),
    activation_energy=19944,
)

experiment_etoh_1 = IdealExperiment(
    name="Romakon-PM102",
    temperature=323.15,
    component=Components.EtOH,
    permeance=Permeance(0.0000282),
    activation_energy=110806,
)

ideal_experiments = IdealExperiments(
    experiments=[
        experiment_h2o_1,
        experiment_etoh_1,
    ]
)

membrane = Membrane(ideal_experiments=ideal_experiments, name="Romakon-PM102")
diffusion_curve = Pervaporation(
    membrane=membrane,
    mixture=Mixtures.H2O_EtOH,
).ideal_diffusion_curve(
    compositions=[
        Composition(p=i / 10, type=CompositionType.weight) for i in range(1, 10)
    ],
    feed_temperature=323.15,
)

measurements = Measurements.from_diffusion_curve_first(diffusion_curve)
best_fit = find_best_curve(measurements, include_zero=False)

membrane = Membrane(
    ideal_experiments=ideal_experiments,
    diffusion_curve_sets=[DiffusionCurveSet("ideal curve", [diffusion_curve])],
    name="Romakon-PM102",
)

conditions = Conditions(
    membrane_area=0.4155,
    initial_feed_temperature=333.15,
    initial_feed_amount=12,
    initial_feed_composition=Composition(p=0.94, type=CompositionType.weight),
)

pervaporation = Pervaporation(membrane, Mixtures.H2O_EtOH)
ideal_model = pervaporation.ideal_isothermal_process(
    conditions=conditions, number_of_steps=1000, delta_hours=0.125
)
non_ideal_model = pervaporation.non_ideal_isothermal_process(
    conditions=conditions,
    diffusion_curves=membrane.diffusion_curve_sets[0],
    number_of_steps=1000,
    delta_hours=0.125,
    n_first=0,
    m_first=0,
    n_second=0,
    m_second=0,
    include_zero=False,
)

x = ideal_model.time
y_ideal = [ideal_model.feed_composition[i].first for i in range(len(x))]
y_non_ideal = [non_ideal_model.feed_composition[i].first for i in range(len(x))]
plt.plot(
    [fc.first for fc in diffusion_curve.feed_compositions],
    [p[0].value for p in diffusion_curve.permeances],
    [fc.first for fc in diffusion_curve.feed_compositions],
    [best_fit(composition.first, diffusion_curve.feed_temperature) for composition in diffusion_curve.feed_compositions])
plt.show()
