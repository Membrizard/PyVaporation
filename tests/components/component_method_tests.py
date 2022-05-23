from component import Component
from utils import AntoineConstants, HeatCapacityConstants

antoine_constants = AntoineConstants(
    a=7.20389,
    b=-1733.926,
    c=-39.485,
)
heat_capacity_constants = HeatCapacityConstants(
    a=32.2,
    b=1.924e-3,
    c=1.055e-5,
    d=-3.596e-9,
)
test_component = Component(
    name="Water",
    molecular_weight=18.02,
    antoine_constants=antoine_constants,
    heat_capacity_constants=heat_capacity_constants,
)


# Covers Antoine Pressure calculation of the components
def test_antoine_pressure():

    # Calculated from https://dx.doi.org/10.1115/1.3687121

    assert abs(test_component.get_antoine_pressure(313) - 7.31934) < 1e-4
    assert abs(test_component.get_antoine_pressure(323) - 12.24821) < 1e-4
    assert abs(test_component.get_antoine_pressure(333) - 19.78961) < 1e-4


# Covers Vaporisation heat calculation from Clapeyron-Clausius equation written with Antoine constants
def test_vaporisation_heat():

    # Molar heat of Vaporisation values for Validation are taken from
    # https://nvlpubs.nist.gov/nistpubs/jres/23/jresv23n2p197_A1b.pdfb

    assert abs(test_component.get_vaporisation_heat(313) - 43.51824) < 3e-1
    assert abs(test_component.get_vaporisation_heat(323) - 42.87114) < 3e-1
    assert abs(test_component.get_vaporisation_heat(333) - 42.43410) < 3e-1


# Covers the calculation of the specific heat with polynomial (n=3) approximation
def test_get_specific_heat():

    # Calculated from http://dx.doi.org/10.1021/je900208n

    assert abs(test_component.get_specific_heat(313) - 33.725516) < 1e-4
    assert abs(test_component.get_specific_heat(323) - 33.800943) < 1e-4
    assert abs(test_component.get_specific_heat(333) - 33.877785) < 1e-4


# Covers the integration of the integral heat required for cooling of the compound,
# including the change in specific heat over a considered temperature range
def test_get_cooling_heat():

    assert abs(test_component.get_cooling_heat(333, 273)-2019.442222) < 2.5e-1
    assert abs(test_component.get_cooling_heat(323, 273)-1681.011314) < 2.5e-1
    assert abs(test_component.get_cooling_heat(313, 273)-1343.342474) < 2.5e-1

