from component import AllComponents
from experiments import IdealExperiments
from mixture import AllMixtures


def test_mixtures_and_components_matching():
    """
    Testing mixtures and components loading.
    Successful loading means all components in mixtures maps on all components in components.
    """
    all_components = AllComponents.load("components.yml")
    all_mixtures = AllMixtures.load("mixtures.yml", all_components=all_components)
    pass


def test_loading_ideal_experiments():
    ideal_experiments = IdealExperiments.from_csv("sample_inputs/IdealExperiment-4.csv")
    pass


def test_antoine_pressure():
    all_components = AllComponents.load("components.yml")
    assert abs(all_components.h2o.get_antoine_pressure(313) - 7.319) < 0.5
