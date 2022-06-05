import numpy
from pytest import fixture

from permeance import Permeance, Units
from component import AllComponents


@fixture
def all_components():
    return AllComponents.load("components.yml")


@fixture
def permeance():
    return Permeance(value=1, units="SI")


def test_converter(permeance, all_components):
    assert (
        permeance.convert(
            to_units=Units().kg_m2_h_kPa, component=all_components.h2o
        ).units
        == "kg/(m2*h*kPa)"
    )
    assert permeance.convert(to_units=Units().GPU).units == "GPU"
    assert permeance.convert(to_units=Units().SI).units == "SI"
    assert permeance.convert(to_units=Units().Barrer).units == "Barrer"
    assert permeance.convert(
        to_units=Units().kg_m2_h_kPa, component=all_components.h2o
    ).value == 1 / (18.02 * 3.6e3)
    assert permeance.convert(to_units=Units().SI).value == 1
    assert permeance.convert(to_units=Units().GPU).value == 1 / 3.35e-10
    assert permeance.convert(to_units=Units().Barrer).value == 1 / 3.35e-16
