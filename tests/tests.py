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
    ideal_experiments = IdealExperiments.from_csv(
        "default_membranes/IdealExperiment-4.csv"
    )
    pass
