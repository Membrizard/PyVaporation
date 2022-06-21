from membrane import Membrane
from mixtures import Mixtures, Composition, CompositionType
from permeance import Permeance, Units
from pervaporation import Pervaporation
from diffusion_curve import DiffusionCurveSet, DiffusionCurve
from experiments import IdealExperiment, IdealExperiments
from components import Components
import matplotlib.pyplot as plt
import numpy
from optimizer import Measurements, find_best_fit
from conditions import Conditions


compositions = [0.99, 0.89, 0.72]

permeance_h2o_50 = [633e-9, 416e-9, 364e-9]
permeance_h2o_60 = [608e-9, 463e-9, 388e-9]
permeance_h2o_70 = [629e-9, 456e-9, 413e-9]
permeance_h2o_80 = [630e-9, 480e-9, 426e-9]

permeance_acetic_acid_50 = [58.4e-9, 3.75e-9, 1.75e-9]
permeance_acetic_acid_60 = [55.4e-9, 3.09e-9, 3.63e-9]
permeance_acetic_acid_70 = [64.7e-9, 5.3e-9, 3.98e-9]
permeance_acetic_acid_80 = [37e-9, 10.3e-9, 4.3e-9]

permeances_h2o = [
    permeance_h2o_50,
    permeance_h2o_60,
    permeance_h2o_70,
    permeance_h2o_80,
]
permeances_acetic_acid = [
    permeance_acetic_acid_50,
    permeance_acetic_acid_60,
    permeance_acetic_acid_70,
    permeance_acetic_acid_80,
]

temperatures = [323.15, 333.15, 343.15, 353.15]
diffusion_curves = []

for t in range(len(temperatures)):
    diffusion_curves.append(
        DiffusionCurve(
            mixture=Mixtures.H2O_AceticAcid,
            membrane_name="Romakon-PM102",
            feed_temperature=temperatures[t],
            feed_compositions=[
                Composition(p=c, type=CompositionType.weight) for c in compositions
            ],
            permeances=[
                (
                    Permeance(value=permeances_h2o[t][i], units=Units.SI),
                    Permeance(value=permeances_acetic_acid[t][i], units=Units.SI),
                )
                for i in range(len(permeances_h2o[t]))
            ],
        )
    )

curve_set = DiffusionCurveSet(
    name_of_the_set="water/acetic acid", diffusion_curves=diffusion_curves
)

measurements_h2o = Measurements.from_diffusion_curves_first(curve_set)
measurements_acetic_acid = Measurements.from_diffusion_curves_second(curve_set)

best_fit_h2o = find_best_fit(data=measurements_h2o)
best_fit_acetic_acid = find_best_fit(data=measurements_acetic_acid)

print(best_fit_h2o)
print(best_fit_acetic_acid)

fig = plt.figure()
ax = fig.add_subplot(projection="3d")

x = [m.x for m in measurements_h2o]
t = [m.t for m in measurements_h2o]
p = [m.p for m in measurements_h2o]
ax.scatter(x, t, p, marker="o")

x_v = numpy.linspace(0.7, 1, num=50)
t_v = numpy.linspace(313.15, 383.15, num=50)
x_fit, t_fit = numpy.meshgrid(x_v, t_v)
p_fit = numpy.array([best_fit_h2o(x_fit[i], t_fit[i]) for i in range(len(x_fit))])
ax.plot_surface(x_fit, t_fit, p_fit, alpha=0.2)

ax.set_xlabel("Water wt.p.")
ax.set_ylabel("Temperature K")
ax.set_zlabel("Water Permeance kg / ( m2 * h * kPa ) ")
fig.suptitle("Best Fit for H2O Illustration", fontsize=10)
plt.show()

fig = plt.figure()
ax = fig.add_subplot(projection="3d")

x = [m.x for m in measurements_acetic_acid]
t = [m.t for m in measurements_acetic_acid]
p = [m.p for m in measurements_acetic_acid]
ax.scatter(x, t, p, marker="o")

x_v = numpy.linspace(0.7, 1, num=50)
t_v = numpy.linspace(313.15, 383.15, num=50)
x_fit, t_fit = numpy.meshgrid(x_v, t_v)
p_fit = numpy.array(
    [best_fit_acetic_acid(x_fit[i], t_fit[i]) for i in range(len(x_fit))]
)
ax.plot_surface(x_fit, t_fit, p_fit, alpha=0.2)

ax.set_xlabel("Water wt.p.")
ax.set_ylabel("Temperature K")
ax.set_zlabel("Water Permeance kg / ( m2 * h * kPa ) ")
fig.suptitle("Best Fit for Acetic Acid Illustration", fontsize=10)
plt.show()


romakon_pm102 = Membrane(name="Romakon_PM102", diffusion_curve_sets=[curve_set])
pervaporation = Pervaporation(membrane=romakon_pm102, mixture=Mixtures.H2O_AceticAcid)

modelled_curve = pervaporation.non_ideal_diffusion_curve(
    diffusion_curve_set=curve_set,
    feed_temperature=373.15,
    initial_feed_composition=Composition(p=0.99, type=CompositionType.weight),
    delta_composition=-0.0054,
    number_of_steps=50,
)

fig = plt.figure()
ax = fig.add_subplot(projection="3d")

x = [m.x for m in measurements_h2o]
t = [m.t for m in measurements_h2o]
p = [m.p for m in measurements_h2o]
ax.scatter(x, t, p, marker="o")

x_v = numpy.linspace(0.7, 1, num=50)
t_v = numpy.linspace(313.15, 383.15, num=50)
x_fit, t_fit = numpy.meshgrid(x_v, t_v)
p_fit = numpy.array([best_fit_h2o(x_fit[i], t_fit[i]) for i in range(len(x_fit))])
ax.plot_surface(x_fit, t_fit, p_fit, alpha=0.2)

x_m = [
    modelled_curve.feed_compositions[i].first
    for i in range(len(modelled_curve.feed_compositions))
]
t_m = [373.15] * (len(x_m))
p_m = [modelled_curve.permeances[i][0].value for i in range(len(x_m))]

ax.scatter(x_m, t_m, p_m, marker="^")

ax.set_xlabel("Water wt.p.")
ax.set_ylabel("Temperature K")
ax.set_zlabel("Water Permeance kg / ( m2 * h * kPa ) ")
fig.suptitle("Modelled Diffusion curve H2O Illustration", fontsize=10)
plt.show()

fig = plt.figure()
ax = fig.add_subplot(projection="3d")

x = [m.x for m in measurements_acetic_acid]
t = [m.t for m in measurements_acetic_acid]
p = [m.p for m in measurements_acetic_acid]
ax.scatter(x, t, p, marker="o")

x_v = numpy.linspace(0.7, 1, num=50)
t_v = numpy.linspace(313.15, 383.15, num=50)
x_fit, t_fit = numpy.meshgrid(x_v, t_v)
p_fit = numpy.array(
    [best_fit_acetic_acid(x_fit[i], t_fit[i]) for i in range(len(x_fit))]
)
ax.plot_surface(x_fit, t_fit, p_fit, alpha=0.2)

x_m = [
    modelled_curve.feed_compositions[i].first
    for i in range(len(modelled_curve.feed_compositions))
]
t_m = [373.15] * len(x_m)
p_m = [modelled_curve.permeances[i][1].value for i in range(len(x_m))]
ax.scatter(x_m, t_m, p_m, marker="^")

ax.set_xlabel("Water wt.p.")
ax.set_ylabel("Temperature K")
ax.set_zlabel("Water Permeance kg / ( m2 * h * kPa ) ")
fig.suptitle("Modelled Diffusion curve Acetic Acid Illustration", fontsize=10)
plt.show()

projected_experiment_conditions = Conditions(
    membrane_area=0.0025,
    initial_feed_temperature=363.15,
    initial_feed_amount=0.050,
    initial_feed_composition=Composition(p=0.90, type=CompositionType.weight),
    permeate_pressure=2,
)

projected_experiment = pervaporation.non_ideal_isothermal_process(
    conditions=projected_experiment_conditions,
    diffusion_curve_set=curve_set,
    number_of_steps=50,
    delta_hours=0.08,
)

composition = [
    projected_experiment.feed_composition[i].second
    for i in range(len(projected_experiment.feed_composition))
]
composition.pop(-1)
plt.plot(projected_experiment.time, composition)
plt.ylabel("Acetic Acid concentration, wt.")
plt.xlabel("Process time, hours")
plt.show()
