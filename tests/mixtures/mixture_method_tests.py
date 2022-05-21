from component import Component
from experiments import IdealExperiments
from mixture import Mixture, Composition, CompositionType, get_nrtl_partial_pressures
from utils import AntoineConstants, HeatCapacityConstants, NRTLParameters

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
test_component_1 = Component(
    name="Water",
    molecular_weight=18.02,
    antoine_constants=antoine_constants,
    heat_capacity_constants=heat_capacity_constants,
)

antoine_constants = AntoineConstants(
    a=7.24677,
    b=-1598.673,
    c=-46.424,
)
heat_capacity_constants = HeatCapacityConstants(
    a=147.815652,
    b=-0.6732612305,
    c=0.001889017424,
    d=0,
)
test_component_2 = Component(
    name="Ethanol",
    molecular_weight=46.07,
    antoine_constants=antoine_constants,
    heat_capacity_constants=heat_capacity_constants,
)

# test_composition_list_molar = [Composition(p=i/10, type=CompositionType.molar) for i in range(10)]
# test_composition_list_weight = [Composition(p=i/10, type=CompositionType.weight) for i in range(10)]

nrtl_params = NRTLParameters(
    g12=-633,
    g21=5823,
    alpha=0.3,
)

test_mixture = Mixture(
    first_component=test_component_1,
    second_component=test_component_2,
    nrtl_params=nrtl_params,
)


def test_composition_to_mol():
    test_composition_1 = Composition(p=0, type=CompositionType("weight"))
    test_composition_2 = Composition(p=1, type=CompositionType("weight"))
    test_composition_3 = Composition(p=0.1, type=CompositionType("weight"))

    assert test_composition_1.to_molar(mixture=test_mixture) == Composition(
        p=0, type=CompositionType("molar")
    )
    assert test_composition_2.to_molar(mixture=test_mixture) == Composition(
        p=1, type=CompositionType("molar")
    )
    assert abs(
        test_composition_3.to_molar(mixture=test_mixture).first - 0.221416
    ) < 1e-4 and test_composition_3.type == CompositionType("molar")


def test_composition_to_weight():
    test_composition_1 = Composition(p=0, type=CompositionType("molar"))
    test_composition_2 = Composition(p=1, type=CompositionType("molar"))
    test_composition_3 = Composition(p=0.221416, type=CompositionType("molar"))
    test_composition_4 = Composition(p=0.630491, type=CompositionType("molar"))

    assert test_composition_1.to_weight(mixture=test_mixture) == Composition(
        p=0, type=CompositionType("weight")
    )
    assert test_composition_2.to_weight(mixture=test_mixture) == Composition(
        p=1, type=CompositionType("weight")
    )
    assert abs(
        test_composition_3.to_weight(mixture=test_mixture).first - 0.1
    ) < 1e-4 and test_composition_3.type == CompositionType("weight")
    assert abs(
        test_composition_4.to_weight(mixture=test_mixture).second - 0.6
    ) < 1e-4 and test_composition_3.type == CompositionType("weight")

