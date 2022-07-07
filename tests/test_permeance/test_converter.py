from pytest import fixture

from pyvaporation.components import Components
from pyvaporation.permeance import Permeance, Units


@fixture
def unit_permeance():
    return Permeance(value=1, units="SI")


@fixture
def real_permeance():
    return Permeance(value=675.64, units="GPU")


def test_converter(unit_permeance):
    assert (
        unit_permeance.convert(
            to_units=Units().kg_m2_h_kPa, component=Components.H2O
        ).units
        == "kg/(m2*h*kPa)"
    )
    assert unit_permeance.convert(to_units=Units().GPU).units == "GPU"
    assert unit_permeance.convert(to_units=Units().SI).units == "SI"
    assert round(
        unit_permeance.convert(
            to_units=Units().kg_m2_h_kPa, component=Components.H2O
        ).value,
        2,
    ) == round(18.02 * 3.6e3, 2)
    assert unit_permeance.convert(to_units=Units().SI).value == 1
    assert unit_permeance.convert(to_units=Units().GPU).value == 1 / 3.35e-10


def test_converter_real_values(real_permeance):
    assert (
        real_permeance.convert(
            to_units=Units().kg_m2_h_kPa, component=Components.H2O
        ).units
        == "kg/(m2*h*kPa)"
    )
    assert real_permeance.convert(to_units=Units().GPU).units == "GPU"
    assert real_permeance.convert(to_units=Units().SI).units == "SI"
    assert (
        abs(
            real_permeance.convert(
                to_units=Units().kg_m2_h_kPa, component=Components.H2O
            ).value
            - 0.0146831
        )
        < 1e-6
    )
    assert (
        abs(real_permeance.convert(to_units=Units().SI).value - 0.00000022634) < 1e-11
    )
    assert real_permeance.convert(to_units=Units().GPU).value == 675.64
