import warnings
from pathlib import Path

import matplotlib
import pytest
from pytest import fixture

from conditions import Conditions
from membrane import Membrane
from mixtures import Composition, CompositionType, Mixtures
from optimizer import Measurements, fit
from permeance import Permeance
from pervaporation import Pervaporation


@fixture
def membrane():

    membrane = Membrane.load(Path("tests/default_membranes/Pervap_2510"))
    return membrane


def auto_close_figures(foo):
    def inner(membrane):
        matplotlib.use(backend="Agg")
        warnings.filterwarnings("ignore")
        foo(membrane)

    return inner


@auto_close_figures
def test_plot_diffusion_curve(membrane):
    diffusion_curve = membrane.diffusion_curve_sets[0][0]
    diffusion_curve.plot(diffusion_curve.partial_fluxes, "Fluxes", curve=False)
    diffusion_curve.plot(diffusion_curve.get_permeances, "Permeances", curve=True)
    diffusion_curve.plot(
        diffusion_curve.get_separation_factor, "Separation Factor", curve=True
    )


@auto_close_figures
def test_plot_pervaporation_function(membrane):
    measurements = Measurements.from_diffusion_curves_first(
        membrane.diffusion_curve_sets[0]
    )

    fit_h2o = fit(measurements)

    fit_h2o.plot()
    fit_h2o.plot(experimental_data=measurements)
    fit_h2o.plot(temperature=333.15)
    fit_h2o.plot(concentration=(0.1, 0.2))
    with pytest.raises(ValueError):
        fit_h2o.plot(experimental_data=measurements, temperature=334.15)
