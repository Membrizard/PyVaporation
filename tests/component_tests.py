from component import AllComponents
from experiments import IdealExperiments
from mixture import AllMixtures

all_components = AllComponents.load("components.yml")


# Covers Antoine Pressure calculation of the components
def test_antoine_pressure():
    assert abs(all_components.h2o.get_antoine_pressure(313) - 7.31934) < 1e-4
    assert abs(all_components.h2o.get_antoine_pressure(323) - 12.24821) < 1e-4
    assert abs(all_components.h2o.get_antoine_pressure(333) - 19.78961) < 1e-4

    assert abs(all_components.meoh.get_antoine_pressure(313) - 35.1986) < 1e-4
    assert abs(all_components.meoh.get_antoine_pressure(323) - 55.2130) < 1e-4
    assert abs(all_components.meoh.get_antoine_pressure(333) - 84.0491) < 1e-4

    assert abs(all_components.etoh.get_antoine_pressure(313) - 17.7708) < 1e-4
    assert abs(all_components.etoh.get_antoine_pressure(323) - 29.2779) < 1e-4
    assert abs(all_components.etoh.get_antoine_pressure(333) - 46.5843) < 1e-4

    assert abs(all_components.ipoh.get_antoine_pressure(313) - 13.7388) < 1e-4
    assert abs(all_components.ipoh.get_antoine_pressure(323) - 23.4142) < 1e-4
    assert abs(all_components.ipoh.get_antoine_pressure(333) - 38.2827) < 1e-4

    assert abs(all_components.dme.get_antoine_pressure(313) - 18.4210) < 1e-4
    assert abs(all_components.dme.get_antoine_pressure(323) - 26.6299) < 1e-4
    assert abs(all_components.dme.get_antoine_pressure(333) - 37.5490) < 1e-4

    assert abs(all_components.dmc.get_antoine_pressure(313) - 14.6761) < 1e-4
    assert abs(all_components.dmc.get_antoine_pressure(323) - 22.9947) < 1e-4
    assert abs(all_components.dmc.get_antoine_pressure(333) - 34.7419) < 1e-4

    assert abs(all_components.mtbe.get_antoine_pressure(313) - 60.1905) < 1e-4
    assert abs(all_components.mtbe.get_antoine_pressure(323) - 85.8647) < 1e-4
    assert abs(all_components.mtbe.get_antoine_pressure(333) - 119.4935) < 1e-4

    assert abs(all_components.etbe.get_antoine_pressure(313) - 31.2979) < 1e-4
    assert abs(all_components.etbe.get_antoine_pressure(323) - 45.9826) < 1e-4
    assert abs(all_components.etbe.get_antoine_pressure(333) - 65.7505) < 1e-4

    assert abs(all_components.cyclohexane.get_antoine_pressure(313) - 24.6260) < 1e-4
    assert abs(all_components.cyclohexane.get_antoine_pressure(323) - 36.2458) < 1e-4
    assert abs(all_components.cyclohexane.get_antoine_pressure(333) - 51.9116) < 1e-4

    assert abs(all_components.benzene.get_antoine_pressure(313) - 24.3752) < 1e-4
    assert abs(all_components.benzene.get_antoine_pressure(323) - 36.1835) < 1e-4
    assert abs(all_components.benzene.get_antoine_pressure(333) - 52.2136) < 1e-4

    assert abs(all_components.toluene.get_antoine_pressure(313) - 7.8949) < 1e-4
    assert abs(all_components.toluene.get_antoine_pressure(323) - 12.2889) < 1e-4
    assert abs(all_components.toluene.get_antoine_pressure(333) - 18.5332) < 1e-4


def test_vaporisation_heat():
    assert abs(all_components.h2o.get_vaporisation_heat(313) - 7.31934) < 1e-4
    assert abs(all_components.h2o.get_vaporisation_heat(323) - 12.24821) < 1e-4
    assert abs(all_components.h2o.get_vaporisation_heat(333) - 19.78961) < 1e-4

    assert abs(all_components.meoh.get_vaporisation_heat(313) - 35.1986) < 1e-4
    assert abs(all_components.meoh.get_vaporisation_heat(323) - 55.2130) < 1e-4
    assert abs(all_components.meoh.get_vaporisation_heat(333) - 84.0491) < 1e-4

    assert abs(all_components.etoh.get_vaporisation_heat(313) - 17.7708) < 1e-4
    assert abs(all_components.etoh.get_vaporisation_heat(323) - 29.2779) < 1e-4
    assert abs(all_components.etoh.get_vaporisation_heat(333) - 46.5843) < 1e-4

    assert abs(all_components.ipoh.get_vaporisation_heat(313) - 13.7388) < 1e-4
    assert abs(all_components.ipoh.get_vaporisation_heat(323) - 23.4142) < 1e-4
    assert abs(all_components.ipoh.get_vaporisation_heat(333) - 38.2827) < 1e-4

    assert abs(all_components.dme.get_vaporisation_heat(313) - 18.4210) < 1e-4
    assert abs(all_components.dme.get_vaporisation_heat(323) - 26.6299) < 1e-4
    assert abs(all_components.dme.get_vaporisation_heat(333) - 37.5490) < 1e-4

    assert abs(all_components.dmc.get_vaporisation_heat(313) - 14.6761) < 1e-4
    assert abs(all_components.dmc.get_vaporisation_heat(323) - 22.9947) < 1e-4
    assert abs(all_components.dmc.get_vaporisation_heat(333) - 34.7419) < 1e-4

    assert abs(all_components.mtbe.get_vaporisation_heat(313) - 60.1905) < 1e-4
    assert abs(all_components.mtbe.get_vaporisation_heat(323) - 85.8647) < 1e-4
    assert abs(all_components.mtbe.get_vaporisation_heat(333) - 119.4935) < 1e-4

    assert abs(all_components.etbe.get_vaporisation_heat(313) - 31.2979) < 1e-4
    assert abs(all_components.etbe.get_vaporisation_heat(323) - 45.9826) < 1e-4
    assert abs(all_components.etbe.get_vaporisation_heat(333) - 65.7505) < 1e-4

    assert abs(all_components.cyclohexane.get_vaporisation_heat(313) - 24.6260) < 1e-4
    assert abs(all_components.cyclohexane.get_vaporisation_heat(323) - 36.2458) < 1e-4
    assert abs(all_components.cyclohexane.get_vaporisation_heat(333) - 51.9116) < 1e-4

    assert abs(all_components.benzene.get_vaporisation_heat(313) - 24.3752) < 1e-4
    assert abs(all_components.benzene.get_vaporisation_heat(323) - 36.1835) < 1e-4
    assert abs(all_components.benzene.get_vaporisation_heat(333) - 52.2136) < 1e-4

    assert abs(all_components.toluene.get_vaporisation_heat(313) - 7.8949) < 1e-4
    assert abs(all_components.toluene.get_vaporisation_heat(323) - 12.2889) < 1e-4
    assert abs(all_components.toluene.get_vaporisation_heat(333) - 18.5332) < 1e-4


def test_get_heat_capacity():
    assert abs(all_components.h2o.get_heat_capacity(300) - 33.629608) < 1e-4


def test_get_specific_heat():
    assert 0 == 0
