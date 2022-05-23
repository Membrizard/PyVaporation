from component import AllComponents
from experiments import IdealExperiments
from mixture import AllMixtures
from membrane import Membrane
from pytest import fixture


@fixture
def all_components():
    return AllComponents.load("components.yml")


@fixture
def component_list(all_components):
    return [all_components.h2o, all_components.meoh, all_components.etoh]


@fixture
def romakon_pm102t():
    ideal_experiments = IdealExperiments.from_csv("default_membranes/RomakonPM102T_Ideal_h2o_alcohol-Table 1.csv")
    return Membrane(ideal_experiments=ideal_experiments, name="Romakon-PM102T")


def test_get_penetrant_data(romakon_pm102t, component_list):
    h2o_experiments = romakon_pm102t.get_penetrant_data(component_list[0])
    meoh_experiments = romakon_pm102t.get_penetrant_data(component_list[1])
    etoh_experiments = romakon_pm102t.get_penetrant_data(component_list[2])

    for i in range(len(h2o_experiments)):
        assert h2o_experiments[i].component == component_list[0]

    for i in range(len(meoh_experiments)):
        assert meoh_experiments[i].component == component_list[1]

    for i in range(len(etoh_experiments)):
        assert etoh_experiments[i].component == component_list[2]


