from pyvaporation.components import Component
from pyvaporation.mixtures import (Composition, CompositionType, Mixture,
                                   get_nrtl_partial_pressures)
from pyvaporation.utils import HeatCapacityConstants, NRTLParameters, VaporPressureConstants

antoine_constants = VaporPressureConstants(
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
    name="h2o",
    molecular_weight=18.02,
    vapour_pressure_constants=antoine_constants,
    heat_capacity_constants=heat_capacity_constants,
)

antoine_constants = VaporPressureConstants(
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
    name="etoh",
    molecular_weight=46.07,
    vapour_pressure_constants=antoine_constants,
    heat_capacity_constants=heat_capacity_constants,
)

test_composition_list_molar = [
    Composition(p=i / 10, type=CompositionType.molar) for i in range(11)
]
test_composition_list_weight = [
    Composition(p=i / 10, type=CompositionType.weight) for i in range(11)
]

nrtl_params = NRTLParameters(
    g12=5823,
    g21=-633,
    alpha12=0.3,
)

test_mixture = Mixture(
    name="H2O_EtOH",
    first_component=test_component_1,
    second_component=test_component_2,
    nrtl_params=nrtl_params,
)


def test_composition_to_mol():
    test_composition_1 = Composition(p=0, type=CompositionType.weight)
    test_composition_2 = Composition(p=1, type=CompositionType.weight)
    test_composition_3 = Composition(p=0.1, type=CompositionType.weight)
    test_composition_4 = Composition(p=0.3, type=CompositionType.weight)

    assert test_composition_1.to_molar(mixture=test_mixture) == Composition(
        p=0, type=CompositionType.molar
    )
    assert test_composition_2.to_molar(mixture=test_mixture) == Composition(
        p=1, type=CompositionType.molar
    )
    assert (
        abs(test_composition_3.to_molar(mixture=test_mixture).first - 0.22122449) < 1e-4
    )
    assert (
        abs(test_composition_4.to_molar(mixture=test_mixture).second - 0.4771704) < 1e-4
    )


def test_composition_to_weight():
    test_composition_1 = Composition(p=0, type=CompositionType.molar)
    test_composition_2 = Composition(p=1, type=CompositionType.molar)
    test_composition_3 = Composition(p=0.22122449, type=CompositionType.molar)
    test_composition_4 = Composition(p=0.63023256, type=CompositionType.molar)

    assert test_composition_1.to_weight(mixture=test_mixture) == Composition(
        p=0, type=CompositionType.weight
    )
    assert test_composition_2.to_weight(mixture=test_mixture) == Composition(
        p=1, type=CompositionType.weight
    )
    assert abs(test_composition_3.to_weight(mixture=test_mixture).first - 0.1) < 1e-4
    assert abs(test_composition_4.to_weight(mixture=test_mixture).second - 0.6) < 1e-4


def test_get_nrtl_partial_pressures_from_molar_composition():
    tested_partial_pressures = [
        get_nrtl_partial_pressures(313, test_mixture, test_composition_list_molar[i])
        for i in range(11)
    ]
    validation_partial_pressures_molar = [
        (0, 17.77081),
        (1.66872, 16.06051),
        (3.06500, 14.49533),
        (4.18580, 13.09037),
        (5.04153, 11.85708),
        (5.65693, 10.80091),
        (6.07206, 9.91318),
        (6.34411, 9.14245),
        (6.55252, 8.28637),
        (6.81318, 6.52717),
        (7.31934, 0),
    ]
    for i in range(11):
        assert (
            abs(
                validation_partial_pressures_molar[i][0]
                - tested_partial_pressures[i][0]
            )
            < 1e-3
        )
    for i in range(11):
        assert (
            abs(
                validation_partial_pressures_molar[i][1]
                - tested_partial_pressures[i][1]
            )
            < 1e-3
        )


def test_get_nrtl_partial_pressures_from_weight_composition():
    tested_partial_pressures = [
        get_nrtl_partial_pressures(313, test_mixture, test_composition_list_weight[i])
        for i in range(11)
    ]
    validation_partial_pressures_weight = [
        (0, 17.77081),
        (3.32577, 14.18329),
        (4.96675, 11.97329),
        (5.76758, 10.58423),
        (6.16626, 9.67241),
        (6.38535, 8.99896),
        (6.53822, 8.35759),
        (6.68342, 7.51991),
        (6.85186, 6.18204),
        (7.06050, 3.89984),
        (7.31934, 0),
    ]
    for i in range(11):
        assert (
            abs(
                validation_partial_pressures_weight[i][0]
                - tested_partial_pressures[i][0]
            )
            < 1e-3
        )
    for i in range(11):
        assert (
            abs(
                validation_partial_pressures_weight[i][1]
                - tested_partial_pressures[i][1]
            )
            < 1e-3
        )
