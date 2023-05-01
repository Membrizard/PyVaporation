from pyvaporation.components import Component
from pyvaporation.mixtures import (
    Composition,
    CompositionType,
    BinaryMixture,
    get_partial_pressures,
)
from pyvaporation.utils import (
    HeatCapacityConstants,
    NRTLParameters,
    VaporPressureConstants,
    UNIQUACConstants,
    UNIQUACParameters,
    UNIQUACBinaryInteractionParameters,
)

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

uniquac_constants = UNIQUACConstants(
    r=0.9200,
    q_geometric=1.4,
    q_interaction=1.00,
)

test_component_1 = Component(
    name="H2O",
    molecular_weight=18.02,
    vapour_pressure_constants=antoine_constants,
    heat_capacity_constants=heat_capacity_constants,
    uniquac_constants=uniquac_constants,
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

uniquac_constants = UNIQUACConstants(
    r=2.10547,
    q_geometric=1.9720,
    q_interaction=0.92,
)

test_component_2 = Component(
    name="EtOH",
    molecular_weight=46.07,
    vapour_pressure_constants=antoine_constants,
    heat_capacity_constants=heat_capacity_constants,
    uniquac_constants=uniquac_constants,
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

uniquac_params = UNIQUACParameters(
    binary_parameters=[
        UNIQUACBinaryInteractionParameters(
            i_component_name="H2O",
            j_component_name="EtOH",
            ij_parameter=(21.127561704493143, -0.9175664931087569),
            ji_parameter=(100.10268878024358, 2.4619377106475753)
        )
    ],
    # alpha_12=21.127561704493143,
    # alpha_21=100.10268878024358,
    # beta_12=-0.9175664931087569,
    # beta_21=2.4619377106475753,
    z=13,
)

test_mixture = BinaryMixture(
    name="H2O_EtOH",
    first_component=test_component_1,
    second_component=test_component_2,
    nrtl_params=nrtl_params,
    uniquac_params=uniquac_params,
).to_mixture()


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
        get_partial_pressures(313, test_mixture, test_composition_list_molar[i])
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
        get_partial_pressures(313, test_mixture, test_composition_list_weight[i])
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


def test_get_uniquac_partial_pressures_from_molar_composition():

    tested_partial_pressures = [
        get_partial_pressures(
            temperature=313,
            mixture=test_mixture,
            composition=test_composition_list_molar[i],
            calculation_type="UNIQUAC",
        )
        for i in range(11)
    ]

    validation_partial_pressures_weight = [
        (0.0, 17.77081102855847),
        (1.7727543232152994, 16.108160331556935),
        (3.07634274576487, 14.65363423233591),
        (4.042497546731469, 13.378776445099259),
        (4.762288693846105, 12.256074257964636),
        (5.302570330576608, 11.253298848587377),
        (5.717207033767349, 10.32067386237128),
        (6.05618949198233, 9.357445263617667),
        (6.375553501678343, 8.113307512353886),
        (6.752795408763098, 5.8539035594786135),
        (7.319336041076023, 0.0),
    ]

    for i in range(11):
        assert (
            abs(
                validation_partial_pressures_weight[i][0]
                - tested_partial_pressures[i][0]
            )
            < 1
        )
    for i in range(11):
        assert (
            abs(
                validation_partial_pressures_weight[i][1]
                - tested_partial_pressures[i][1]
            )
            < 1
        )


def test_get_uniquac_partial_pressures_from_weight_composition():

    tested_partial_pressures = [
        get_partial_pressures(
            temperature=313,
            mixture=test_mixture,
            composition=test_composition_list_weight[i],
            calculation_type="UNIQUAC",
        )
        for i in range(11)
    ]

    validation_partial_pressures_weight = [
        (0.0, 17.77081102855847),
        (3.3060551887742338, 14.368949786545734),
        (4.698897748463025, 12.363107937595991),
        (5.4064961856347935, 11.036649571998714),
        (5.82532411896513, 10.039774236003971),
        (6.115906467217008, 9.156764434038331),
        (6.352970243831038, 8.217415786032971),
        (6.574748311471986, 7.053789671555542),
        (6.8022687199985175, 5.472027024051633),
        (7.047875750679181, 3.227809119819075),
        (7.319336041076023, 0.0),
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
