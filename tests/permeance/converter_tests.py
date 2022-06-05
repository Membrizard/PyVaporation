import numpy
from pytest import fixture

from permeance import Permeance, Units
from component import AllComponents


@fixture
def all_components():
    return AllComponents.load("components.yml")


@fixture
def unit_permeance():
    return Permeance(value=1, units="SI")


@fixture
def real_permeance():
    return Permeance(value=675.64, units="GPU")


def test_converter(unit_permeance, all_components):
    assert (
            unit_permeance.convert(
                to_units=Units().kg_m2_h_kPa, component=all_components.h2o
            ).units
            == "kg/(m2*h*kPa)"
    )
    assert unit_permeance.convert(to_units=Units().GPU).units == "GPU"
    assert unit_permeance.convert(to_units=Units().SI).units == "SI"
    assert round(unit_permeance.convert(
        to_units=Units().kg_m2_h_kPa, component=all_components.h2o
    ).value, 2) == round(18.02 * 3.6e3, 2)
    assert unit_permeance.convert(to_units=Units().SI).value == 1
    assert unit_permeance.convert(to_units=Units().GPU).value == 1 / 3.35e-10


def test_converter_real_values(real_permeance, all_components):
    assert (
        real_permeance.convert(
            to_units=Units().kg_m2_h_kPa, component=all_components.h2o
        ).units
        == "kg/(m2*h*kPa)"
    )
    assert real_permeance.convert(to_units=Units().GPU).units == "GPU"
    assert real_permeance.convert(to_units=Units().SI).units == "SI"
    assert abs(real_permeance.convert(
        to_units=Units().kg_m2_h_kPa, component=all_components.h2o
    ).value - 0.0146831) < 1e-6
    assert abs(real_permeance.convert(to_units=Units().SI).value - 0.00000022634) < 1e-11
    assert real_permeance.convert(to_units=Units().GPU).value == 675.64
