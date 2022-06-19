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


pervaporation = Pervaporation(pervap_4100, Mixtures.H2O_EtOH)
modelled_curve = pervaporation.non_ideal_diffusion_curve(
    pervap_4100.diffusion_curve_sets[0],
    368.15,
    experimental_50.feed_compositions,
)
print(
    [
        modelled_curve.permeances[i][1].value
        for i in range(len(modelled_curve.permeances))
    ]
)
x = compositions_50
y_ideal_EtOH = [experimental_50.permeances[i][1].value for i in range(len(x))]
y_non_ideal_EtOH = [modelled_curve.permeances[i][1].value for i in range(len(x))]
y_ideal_H2O = [experimental_50.permeances[i][0].value for i in range(len(x))]
y_non_ideal_H2O = [modelled_curve.permeances[i][0].value for i in range(len(x))]
plt.plot(x, y_ideal_EtOH, x, y_non_ideal_EtOH, x, y_ideal_H2O, x, y_non_ideal_H2O)
plt.legend(["Experimental EtOH_50", "Modelled EtOH_50", "Experimental H2O_50", "Modelled H2O_50"])
plt.show()

compositions_25 = [
        0.005373433856,
        0.03086565059,
        0.03764549547,
        0.1789374627,
        0.2239556327,
        0.2570412757,
    ]

permeances_25_H2O = [
        1212.734594,
        1561.012113,
        1703.565186,
        3202.285314,
        3329.107307,
        3598.01802380722,
    ]

permeances_25_EtOH = [
        1.778198976,
        2.426177206,
        2.522262532,
        10.21010056,
        11.92617288,
        13.01528204,
    ]

compositions_25.reverse()
permeances_25_H2O.reverse()
permeances_25_EtOH.reverse()

permeances_25_H2O = [Permeance(permeance,Units.GPU).convert(Units.kg_m2_h_kPa,Components.H2O).value for permeance in permeances_25_H2O]
permeances_25_EtOH = [Permeance(permeance,Units.GPU).convert(Units.kg_m2_h_kPa,Components.EtOH).value for permeance in permeances_25_EtOH]


modelled_curve_25 = pervaporation.non_ideal_diffusion_curve(
        diffusion_curves=pervaporation.membrane.diffusion_curve_sets[0],
        feed_temperature=368.15,
        compositions=[Composition(p=composition, type=CompositionType.weight) for composition in compositions_25],
        initial_permeances=(Permeance(value=3598.01802380722, units=Units.GPU),
                             Permeance(value=13.01528204, units=Units.GPU)),
    )

x = compositions_25
y_ideal_EtOH = [permeances_25_EtOH[i] for i in range(len(x))]
y_non_ideal_EtOH = [modelled_curve_25.permeances[i][1].value for i in range(len(x))]
y_ideal_H2O = [permeances_25_H2O[i] for i in range(len(x))]
y_non_ideal_H2O = [modelled_curve_25.permeances[i][0].value for i in range(len(x))]
plt.plot(x, y_ideal_EtOH, x, y_non_ideal_EtOH, x, y_ideal_H2O, x, y_non_ideal_H2O)
plt.legend(["Experimental EtOH_25", "Modelled EtOH_25", "Experimental H2O_25", "Modelled H2O_25"])
plt.show()