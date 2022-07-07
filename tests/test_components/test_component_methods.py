from pyvaporation.components import Component
from pyvaporation.utils import HeatCapacityConstants, VaporPressureConstants

vapor_pressure_constants_antoine = VaporPressureConstants(
    a=7.20389,
    b=-1733.926,
    c=-39.485,
)

vapor_pressure_constants_frost = VaporPressureConstants(
    a=16.5191,
    b=-3937.6553,
    c=-190231.9062,
    type="frost",
)

heat_capacity_constants = HeatCapacityConstants(
    a=32.2,
    b=1.924e-3,
    c=1.055e-5,
    d=-3.596e-9,
)
test_component = Component(
    name="H2O_antoine",
    molecular_weight=18.02,
    vapour_pressure_constants=vapor_pressure_constants_antoine,
    heat_capacity_constants=heat_capacity_constants,
)

test_component_2 = Component(
    name="H2O_frost",
    molecular_weight=18.02,
    vapour_pressure_constants=vapor_pressure_constants_frost,
    heat_capacity_constants=heat_capacity_constants,
)


# Covers Antoine Pressure calculation of the test_components
def test_antoine_pressure():

    # Calculated from https://dx.doi.org/10.1115/1.3687121

    assert abs(test_component.get_vapor_pressure(313) - 7.31934) < 1e-4
    assert abs(test_component.get_vapor_pressure(323) - 12.24821) < 1e-4
    assert abs(test_component.get_vapor_pressure(333) - 19.78961) < 1e-4


def test_frost_pressure():
    # Pressures and temperatures for validation are taken from http://www.ddbst.com/en/EED/PCP/VAP_C174.php
    validation_temperatures = [293.15, 313.15, 323.15, 333.15, 343.15, 412.75]
    validation_pressures = [2.400, 7.333, 12.332, 19.932, 31.224, 357.0]
    for i in range(len(validation_pressures)):
        assert (
            abs(
                test_component_2.get_vapor_pressure(validation_temperatures[i])
                - validation_pressures[i]
            )
            < validation_pressures[i] * 0.02
        )
        assert (
            abs(
                test_component_2.get_vapor_pressure(validation_temperatures[i])
                - test_component.get_vapor_pressure(validation_temperatures[i])
            )
            < abs(
                test_component_2.get_vapor_pressure(validation_temperatures[i])
                + test_component.get_vapor_pressure(validation_temperatures[i])
            )
            * 0.03
            / 2
        )


def test_frost_vaporisation_heat():
    validation_temperatures = [298.15, 373.10, 373.15, 373.20]
    validation_enthalpies = [43.950, 40.7, 40.599, 41.110]

    for i in range(len(validation_temperatures)):
        assert (
            abs(
                test_component_2.get_vaporisation_heat(validation_temperatures[i])
                - validation_enthalpies[i]
            )
            < validation_enthalpies[i] * 0.02
        )


def test_antoine_vaporisation_heat():
    validation_temperatures = [298.15, 373.10, 373.15, 373.20]
    validation_enthalpies = [43.950, 40.7, 40.599, 41.110]

    for i in range(len(validation_temperatures)):
        assert (
            abs(
                test_component.get_vaporisation_heat(validation_temperatures[i])
                - validation_enthalpies[i]
            )
            < validation_enthalpies[i] * 0.03
        )


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

    assert abs(test_component.get_cooling_heat(333, 273) - 2019.442222) < 2.5e-1
    assert abs(test_component.get_cooling_heat(323, 273) - 1681.011314) < 2.5e-1
    assert abs(test_component.get_cooling_heat(313, 273) - 1343.342474) < 2.5e-1
