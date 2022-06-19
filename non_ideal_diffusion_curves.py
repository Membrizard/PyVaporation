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

pervap_4100 = Membrane(
    name=experimental_50.membrane_name,
    diffusion_curve_sets=[
        DiffusionCurveSet(
            name_of_the_set="Various Initial Feed", diffusion_curves=[experimental_50]
        )
    ],
)

compositions_6 = [0.01703476704,
                  0.0194755112,
                  0.04876444107,
                  0.05473070456,
                  0.05961219287]

permeances_6_H2O =[
746.3089286,
838.5364554,
1038.232656,
1122.096546,
1133.045025
]

permeances_6_H2O.reverse()
compositions_6.reverse()
compositions_6 = [Composition(composition, CompositionType.weight) for composition in compositions_6]

pervaporation = Pervaporation(pervap_4100, Mixtures.H2O_EtOH)

modelled_curve = pervaporation.non_ideal_diffusion_curve(diffusion_curves=pervaporation.membrane.diffusion_curve_sets[0],
                                                         feed_temperature=368.15,
                                                         compositions=compositions_6,
                                                         initial_permeances=(Permeance(1133.045025, Units.GPU),
                                                                             Permeance(0.2100162934, Units.GPU)))

x = [composition.first for composition in compositions_6]
y_ideal = [Permeance(permeance, Units.GPU).convert(Units.kg_m2_h_kPa, Components.H2O).value for permeance in permeances_6_H2O]
y_non_ideal = [modelled_curve.permeances[i][0].value for i in range(len(x))]
plt.plot(x, y_ideal, x, y_non_ideal)
plt.legend(["experimental", "modelled"])
plt.show()
