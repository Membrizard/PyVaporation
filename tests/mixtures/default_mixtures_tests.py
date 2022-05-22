from component import Component, AllComponents
from mixture import (
    Mixture,
    Composition,
    CompositionType,
    AllMixtures,
    get_nrtl_partial_pressures,
)
from utils import AntoineConstants, HeatCapacityConstants, NRTLParameters

all_components = AllComponents.load("components.yml")
all_mixtures = AllMixtures.load("mixtures.yml", all_components)


def test_constants_h2o_etoh():
    # Experimental data for validation is taken from https://doi.org/10.1021/je00019a033
    # NRTL Constants are taken from Identification of Best Model for Equilibrium Data of Ethanol-Water Mixture
    # Bilel Hadrich and Nabil Kechaou
    # June 2010, Volume 4, No.6 (Serial No.31)Journal of Chemistry and Chemical Engineering, ISSN 1934-7375, USA
    test_mixture = Mixture.from_dict(d="h2o_etoh", all_components=all_components)
    validation_compositions = [
        Composition(p=0.82440, type=CompositionType("molar")),
        Composition(p=0.62270, type=CompositionType("molar")),
        Composition(p=0.40308, type=CompositionType("molar")),
        Composition(p=0.09690, type=CompositionType("molar")),
    ]

    validation_pressures = [
        (10.9512364, 12.7117636),
        (9.9939294, 16.4870706),
        (8.5117929, 19.5892071),
        (2.7611485, 26.7698515),
    ]
    tested_partial_pressures = [
        get_nrtl_partial_pressures(
            temperature=323.15,
            mixture=test_mixture,
            composition=validation_compositions[i],
        )
        for i in range(5)
    ]

    for i in range(5):
        assert abs(tested_partial_pressures[i][0] - validation_pressures[i][0]) < 1e-2
        assert abs(tested_partial_pressures[i][1] - validation_pressures[i][1]) < 1e-2


def test_constants_h2o_ipoh():
    # Experimental data for validation is taken from http://www.ddbst.com/en/EED/VLE/VLE%202-Propanol%3BWater.php
    # Dunlop J.G.: Vapor-Liquid Equilibrium Data. Master's Thesis (1948)
    # NRTL Constants are taken from https://doi.org/10.1021/je960108n
    test_mixture = Mixture.from_dict(d="h2o_ipoh", all_components=all_components)
    validation_compositions = [
        Composition(p=0.9800, type=CompositionType("molar")),
        Composition(p=0.7610, type=CompositionType("molar")),
        Composition(p=0.3030, type=CompositionType("molar")),
        Composition(p=0.06800, type=CompositionType("molar")),
    ]

    validation_temperature_list = [363.95, 354.25, 353.25, 354.35]

    validation_pressures = [
        (77.92277, 23.407230),
        (44.99052, 56.339480),
        (31.10831, 70.221690),
        (10.03167, 91.298330),
    ]
    tested_partial_pressures = [
        get_nrtl_partial_pressures(
            temperature=validation_temperature_list[i],
            mixture=test_mixture,
            composition=validation_compositions[i],
        )
        for i in range(5)
    ]

    for i in range(5):
        assert abs(tested_partial_pressures[i][0] - validation_pressures[i][0]) < 1e-2
        assert abs(tested_partial_pressures[i][1] - validation_pressures[i][1]) < 1e-2


def test_constants_etoh_etbe():
    # Experimental data for validation and NRTL constants are taken from
    # Isothermal vapor-liquid equilibria for binary and ternary systems containing ethyl tert-butyl ether,
    # ethanol, benzene, and toluene at 313.15 K
    # Oh, JH; Park, SJ
    # Journal of Industrial and Engineering Chemistry, 2005

    test_mixture = Mixture.from_dict(d="etoh_etbe", all_components=all_components)
    validation_compositions = [
        Composition(p=0.9007, type=CompositionType("molar")),
        Composition(p=0.5026, type=CompositionType("molar")),
        Composition(p=0.1994, type=CompositionType("molar")),
        Composition(p=0.0204, type=CompositionType("molar")),
    ]

    validation_pressures = [
        (16.29117, 7.64883),
        (11.612918, 21.397082),
        (7.273212, 27.066788),
        (1.128097, 31.381903),
    ]
    tested_partial_pressures = [
        get_nrtl_partial_pressures(
            temperature=313.15,
            mixture=test_mixture,
            composition=validation_compositions[i],
        )
        for i in range(5)
    ]

    for i in range(5):
        assert abs(tested_partial_pressures[i][0] - validation_pressures[i][0]) < 1e-2
        assert abs(tested_partial_pressures[i][1] - validation_pressures[i][1]) < 1e-2


def test_constants_meoh_toluene():
    # Experimental data for validation is taken from https://doi.org/10.1021/je60021a018
    # NRTL constants are taken from Understanding Distillation Using Column Profile Maps, F
    # irst Edition. Daniel Beneke, Mark Peters, David Glasser, and Diane Hildebrandt.
    # 2013 by John Wiley & Sons, Inc. Published 2013 by John Wiley & Sons, Inc
    # https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781118477304.app2

    test_mixture = Mixture.from_dict(d="meoh_toluene", all_components=all_components)
    validation_compositions = [
        Composition(p=0.1320, type=CompositionType("molar")),
        Composition(p=0.4390, type=CompositionType("molar")),
        Composition(p=0.8300, type=CompositionType("molar")),
        Composition(p=0.9740, type=CompositionType("molar")),
    ]

    validation_pressures = [
        (81.161325,	20.163675),
        (83.8971,	17.4279),
        (87.74745,	13.57755),
        (96.968025,	4.356975),
    ]

    validation_temperature_list = [342.70, 338.10, 336.70, 337.10]

    tested_partial_pressures = [
        get_nrtl_partial_pressures(
            temperature=validation_temperature_list[i],
            mixture=test_mixture,
            composition=validation_compositions[i],
        )
        for i in range(5)
    ]

    for i in range(5):
        assert abs(tested_partial_pressures[i][0] - validation_pressures[i][0]) < 1e-2
        assert abs(tested_partial_pressures[i][1] - validation_pressures[i][1]) < 1e-2


def test_constants_meoh_mtbe():
    # Experimental data for validation and NRTL Parameters are taken from
    # https://doi.org/10.1002/1521-4125(20020709)25:7<729::AID-CEAT729>3.0.CO;2-B

    test_mixture = Mixture.from_dict(d="meoh_mtbe", all_components=all_components)
    validation_compositions = [
        Composition(p=0.9007, type=CompositionType("molar")),
        Composition(p=0.5026, type=CompositionType("molar")),
        Composition(p=0.1994, type=CompositionType("molar")),
        Composition(p=0.0204, type=CompositionType("molar")),
    ]

    validation_pressures = [
        (16.29117, 7.64883),
        (11.612918, 21.397082),
        (7.273212, 27.066788),
        (1.128097, 31.381903),
    ]
    tested_partial_pressures = [
        get_nrtl_partial_pressures(
            temperature=313.15,
            mixture=test_mixture,
            composition=validation_compositions[i],
        )
        for i in range(5)
    ]

    for i in range(5):
        assert abs(tested_partial_pressures[i][0] - validation_pressures[i][0]) < 1e-2
        assert abs(tested_partial_pressures[i][1] - validation_pressures[i][1]) < 1e-2


