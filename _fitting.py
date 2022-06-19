from optimizer import Measurements, fit, find_best_fit
import matplotlib.pyplot as plt
from components import Components
from conditions import Conditions
from experiments import IdealExperiment, IdealExperiments
from membrane import Membrane
from mixtures import Mixtures, Composition, CompositionType, get_nrtl_partial_pressures
from permeance import Permeance, Units
from pervaporation import Pervaporation
from diffusion_curve import DiffusionCurveSet, DiffusionCurve

compositions_50 = [
    0.02408580572,
    0.2326338341,
    0.2532445626,
    0.2754824538,
    0.4110793513,
    0.4284357542,
    0.4517584206,
    0.4745386994,
]

permeances_50_H2O = [
    2688.778689,
    5623.966085,
    5790.19941,
    6078.245449,
    7525.770499,
    7823.818024,
    7525.770499,
    7525.770499,
]
permeances_50_EtOH = [
    23.76315076,
    107.036062,
    120.2633877,
    132.5265399,
    207.1487026,
    215.3525347,
    205.1470494,
    205.1470494,
]
experimental_50 = DiffusionCurve(
    mixture=Mixtures.H2O_EtOH,
    membrane_name="Pervap 4100",
    feed_temperature=368.15,
    feed_compositions=[
        Composition(p=composition, type=CompositionType.weight)
        for composition in compositions_50
    ],
    permeate_pressure=1,
    permeances=[
        (
            Permeance(value=permeances_50_H2O[i], units=Units.GPU),
            Permeance(value=permeances_50_EtOH[i], units=Units.GPU),
        )
        for i in range(len(compositions_50))
    ],
)
measurements_etoh = Measurements.from_diffusion_curve_second(experimental_50)

fit_etoh = find_best_fit(data=measurements_etoh, n=0, m=1)
print(fit_etoh)
x = compositions_50
y_ideal = [
    Permeance(permeance, Units.GPU).convert(Units.kg_m2_h_kPa, Components.EtOH).value
    for permeance in permeances_50_EtOH
]
y_non_ideal = [fit_etoh(composition, 368.15) for composition in compositions_50]
plt.plot(x, y_ideal, x, y_non_ideal)
plt.legend(["experiment", "Model"])
plt.show()
