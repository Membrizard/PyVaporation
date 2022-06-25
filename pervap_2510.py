from mixtures import Mixtures, Composition, CompositionType, get_nrtl_partial_pressures
from components import Components
from diffusion_curve import DiffusionCurve, DiffusionCurveSet
from membrane import Membrane
from mixtures import Composition, CompositionType, Mixtures
from permeance import Permeance
from optimizer import Measurements, find_best_fit
from pervaporation import Pervaporation
from conditions import Conditions
import matplotlib.pyplot as plt
import numpy

compositions_mol_100 = [
    0.0909,
    0.1354,
    0.2081,
    0.2604,
    0.3164,
]
pervap_2510_mol_flux_water_100 = [29.7928, 46.5392, 95.3431, 128.3574, 188.1660]
pervap_2510_mol_flux_ipoh_100 = [0.4597, 0.6290, 0.8226, 1.0060, 1.3831]

compositions_mol_90 = [
    0.1251,
    0.1825,
    0.2379,
    0.2939,
    0.3438,
]
pervap_2510_mol_flux_water_90 = [30.2713, 46.5392, 76.4435, 106.8263, 134.5775]
pervap_2510_mol_flux_ipoh_90 = [0.25, 0.3266, 0.5443, 0.7238, 0.9516]

compositions_mol_80 = [0.201, 0.2461, 0.3007, 0.3305, 0.3643]
pervap_2510_mol_flux_water_80 = [36.49139, 57.5440, 68.0703, 81.9459, 103.4770]
pervap_2510_mol_flux_ipoh_80 = [0.1734, 0.2823, 0.3468, 0.3790, 0.4516]

compositions_mol_70 = [0.2604, 0.2966, 0.3383, 0.3506]
pervap_2510_mol_flux_water_70 = [34.5775, 42.7115, 48.6923, 80.5105]
pervap_2510_mol_flux_ipoh_70 = [0.0968, 0.1169, 0.1976, 0.1814]

compositions_mol_60 = [0.2932, 0.3349, 0.3414, 0.3691]
pervap_2510_mol_flux_water_60 = [29.3143, 32.6636, 36.9699, 47.0177]
pervap_2510_mol_flux_ipoh_60 = [0.0524, 0.0685, 0.0847, 0.0726]

all_compositions = [
    compositions_mol_60,
    compositions_mol_70,
    compositions_mol_80,
    compositions_mol_90,
    compositions_mol_100,
]

all_fluxes_h2o = [
    pervap_2510_mol_flux_water_60,
    pervap_2510_mol_flux_water_70,
    pervap_2510_mol_flux_water_80,
    pervap_2510_mol_flux_water_90,
    pervap_2510_mol_flux_water_100,
]

all_fluxes_ipoh = [
    pervap_2510_mol_flux_ipoh_60,
    pervap_2510_mol_flux_ipoh_70,
    pervap_2510_mol_flux_ipoh_80,
    pervap_2510_mol_flux_ipoh_90,
    pervap_2510_mol_flux_ipoh_100,
]

temperatures = [333.15, 343.15, 353.15, 363.15, 373.15]
curves = []

for i in range(len(temperatures)):
    curves.append(
        DiffusionCurve(
            mixture=Mixtures.H2O_iPOH,
            membrane_name="Pervap 2510",
            feed_temperature=temperatures[i],
            feed_compositions=[
                Composition(p=c, type=CompositionType.molar)
                for c in all_compositions[i]
            ],
            partial_fluxes=[
                (
                    all_fluxes_h2o[i][j] * Components.H2O.molecular_weight / 1000,
                    all_fluxes_ipoh[i][j] * Components.iPOH.molecular_weight / 1000,
                )
                for j in range(len(all_compositions[i]))
            ],
        )
    )

curve_set = DiffusionCurveSet(
    name="water - isopropanol 5-30 mol% ", diffusion_curves=curves
)

pervap_2510 = Membrane(
    name="Pervap 2510",
    diffusion_curve_sets=[curve_set],
)

pervaporation = Pervaporation(membrane=pervap_2510, mixture=Mixtures.H2O_iPOH)

conditions = Conditions(
        membrane_area=30,
        initial_feed_temperature=373.15,
        initial_feed_amount=1000,
        initial_feed_composition=Composition(p=0.15, type=CompositionType.weight),
    )

modelled_process = pervaporation.non_ideal_non_isothermal_process(
        conditions=conditions,
        diffusion_curve_set=pervap_2510.diffusion_curve_sets[0],
        number_of_steps=10,
        delta_hours=0.1,
        initial_permeances=(Permeance(0.056), Permeance(0.00014)),
    )

module_area_weight_fraction = [1.830, 5.297, 9.718, 14.660, 19.298, 23.763]
validation_water_wt_fraction = [0.1441, 0.1334, 0.1227, 0.1125, 0.1048, 0.0988]

module_area_temperature = [1.645, 3.575, 6.776, 10.592, 14.803, 19.978, 27.566]
validation_temperature = [98.679, 96.710, 94.201, 91.548, 89.061, 86.408, 83.312]

module_area_flux = [0.87, 3.78, 7.46, 12.12, 18.88, 23.71, 28.41]
validation_flux = [4.1880, 3.4200, 2.7480, 2.1360, 1.6080, 1.3440, 1.1280]

pressure = get_nrtl_partial_pressures(
    373.15, Mixtures.H2O_iPOH, Composition(p=0.15, type=CompositionType.weight)
)
print(validation_flux[0] / pressure[0])

plt.plot([c*30 for c in modelled_process.time], [c.first for c in modelled_process.feed_composition])
plt.plot(module_area_weight_fraction, validation_water_wt_fraction, label="literature")
plt.legend(["model", "experiment"])
plt.show()

plt.plot(module_area_temperature, [t + 273.15 for t in validation_temperature])
plt.plot([c*30 for c in modelled_process.time], modelled_process.feed_temperature)
plt.legend(["experiment", "model"])
plt.show()

plt.plot(module_area_flux, validation_flux)
plt.plot([c*30 for c in modelled_process.time], [fluxes[0] for fluxes in modelled_process.partial_fluxes])
plt.legend(["experiment", "model"])
plt.show()

measurements_h2o = Measurements.from_diffusion_curves_first(curve_set)
measurements_ipoh = Measurements.from_diffusion_curves_second(curve_set)

fit_h2o = find_best_fit(measurements_h2o)
fit_ipoh = find_best_fit(measurements_ipoh)

print(fit_h2o)
print(fit_ipoh)

fig = plt.figure()
ax = fig.add_subplot(projection="3d")

x = [m.x for m in measurements_h2o]
t = [m.t for m in measurements_h2o]
p = [m.p for m in measurements_h2o]
ax.scatter(x, t, p, marker="o")

x_v = numpy.linspace(0, 0.4, num=50)
t_v = numpy.linspace(323.15, 383.15, num=50)
x_fit, t_fit = numpy.meshgrid(x_v, t_v)
p_fit = numpy.array([fit_h2o(x_fit[i], t_fit[i]) for i in range(len(x_fit))])
ax.plot_surface(x_fit, t_fit, p_fit, alpha=0.2)

ax.set_xlabel("Water wt.p.")
ax.set_ylabel("Temperature K")
ax.set_zlabel("Water Permeance kg / ( m2 * h * kPa ) ")
fig.suptitle("Best Fit for water Illustration", fontsize=10)
plt.show()

fig = plt.figure()
ax = fig.add_subplot(projection="3d")

x = [m.x for m in measurements_ipoh]
t = [m.t for m in measurements_ipoh]
p = [m.p for m in measurements_ipoh]
ax.scatter(x, t, p, marker="o")

x_v = numpy.linspace(0, 0.4, num=50)
t_v = numpy.linspace(323.15, 383.15, num=50)
x_fit, t_fit = numpy.meshgrid(x_v, t_v)
p_fit = numpy.array([fit_ipoh(x_fit[i], t_fit[i]) for i in range(len(x_fit))])
ax.plot_surface(x_fit, t_fit, p_fit, alpha=0.2)

ax.set_xlabel("Water wt.p.")
ax.set_ylabel("Temperature K")
ax.set_zlabel("iPOH Permeance kg / ( m2 * h * kPa ) ")
fig.suptitle("Best Fit for isopropanol Illustration", fontsize=10)
plt.show()
