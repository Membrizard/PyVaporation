from membrane import Membrane
from mixtures import Mixtures, Composition, CompositionType
from permeance import Permeance, Units
from pervaporation import Pervaporation
from diffusion_curve import DiffusionCurveSet, DiffusionCurve
import matplotlib.pyplot as plt
import numpy
from experiments import IdealExperiment, IdealExperiments
from components import Components
from optimizer import Measurements, find_best_fit

compositions = [1, 0.883, 0.756, 0.658, 0.490, 0.430, 0.401, 0.280, 0.240, 0.044]

flux_h2o_40 = [
        0.7848,
        0.3480,
        0.2566,
        0.4406,
        0.3340,
        0.2773,
        0.2816,
        0.2292,
        0.1066,
        0.0247,
    ]

flux_etoh_40 = [
        0,
        0.0017,
        0.0013,
        0.0136,
        0.0684,
        0.0568,
        0.0577,
        0.0573,
        0.0457,
        0.0165,
    ]

spi_curve = DiffusionCurve(
        mixture=Mixtures.H2O_EtOH,
        membrane_name="SPI 255 dense",
        feed_temperature=313.15,
        feed_compositions=[
            Composition(p=c, type=CompositionType.weight) for c in compositions
        ],
        partial_fluxes=[
            (flux_h2o_40[i], flux_etoh_40[i]) for i in range(len(compositions))
        ],
    )

measurements_h2o = Measurements.from_diffusion_curve_first(spi_curve)
measurements_etoh = Measurements.from_diffusion_curve_second(spi_curve)

fit_h2o = find_best_fit(measurements_h2o, n=9)
fit_etoh = find_best_fit(measurements_etoh)

print(fit_h2o)

p_h2o = [permeance[0].value for permeance in spi_curve.permeances]
x_m = numpy.arange(0, 1.05, 0.05)
plt.plot(compositions, p_h2o)
plt.plot(x_m, [fit_h2o(x, 313.15) for x in x_m])
plt.show()