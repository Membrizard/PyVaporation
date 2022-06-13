import numpy
from pytest import fixture

from components import Components
from mixtures import (
    Mixtures,
    Composition,
    CompositionType,
    get_nrtl_partial_pressures,
)


def test_constants_h2o_meoh():
    # Experimental data for validation is taken from https://doi.org/10.1021/je00019a033
    # NRTL Constants are taken from https://doi.org/10.1002/9781118477304.app2
    test_mixture = Mixtures.H2O_MeOH
    validation_compositions = [
        Composition(p=0.7530, type=CompositionType("molar")),
        Composition(p=0.5972, type=CompositionType("molar")),
        Composition(p=0.3855, type=CompositionType("molar")),
        Composition(p=0.2270, type=CompositionType("molar")),
    ]

    validation_pressures = [
        (9.580151, 19.538849),
        (7.8695188, 27.4514812),
        (6.0172119, 36.0317881),
        (3.976056, 43.357944),
    ]
    tested_partial_pressures = [
        get_nrtl_partial_pressures(
            temperature=323.15,
            mixture=test_mixture,
            composition=validation_compositions[i],
        )
        for i in range(4)
    ]

    rmsd_1 = 0
    rmsd_2 = 0

    for i in range(4):
        rmsd_1 = (
            tested_partial_pressures[i][0] - validation_pressures[i][0]
        ) ** 2 / validation_pressures[i][0] ** 2 + rmsd_1

        rmsd_2 = (
            tested_partial_pressures[i][1] - validation_pressures[i][1]
        ) ** 2 / validation_pressures[i][1] ** 2 + rmsd_2
        assert (
            abs(tested_partial_pressures[i][0] - validation_pressures[i][0])
            < validation_pressures[i][0] * 0.05
        )
        assert (
            abs(tested_partial_pressures[i][1] - validation_pressures[i][1])
            < validation_pressures[i][1] * 0.05
        )
    assert numpy.sqrt(rmsd_1 / 4) < 0.033
    assert numpy.sqrt(rmsd_2 / 4) < 0.033


def test_constants_h2o_etoh():
    # Experimental data for validation is taken from https://doi.org/10.1021/je00019a033
    # NRTL Constants are taken from Identification of Best Model for Equilibrium Data of Ethanol-Water Mixture
    # Bilel Hadrich and Nabil Kechaou
    # June 2010, Volume 4, No.6 (Serial No.31)Journal of Chemistry and Chemical Engineering, ISSN 1934-7375, USA
    test_mixture = Mixtures.H2O_EtOH
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
        for i in range(4)
    ]

    rmsd_1 = 0
    rmsd_2 = 0

    for i in range(4):
        rmsd_1 = (
            tested_partial_pressures[i][0] - validation_pressures[i][0]
        ) ** 2 / validation_pressures[i][0] ** 2 + rmsd_1

        rmsd_2 = (
            tested_partial_pressures[i][1] - validation_pressures[i][1]
        ) ** 2 / validation_pressures[i][1] ** 2 + rmsd_2
        assert (
            abs(tested_partial_pressures[i][0] - validation_pressures[i][0])
            < validation_pressures[i][0] * 0.05
        )
        assert (
            abs(tested_partial_pressures[i][1] - validation_pressures[i][1])
            < validation_pressures[i][1] * 0.05
        )
    assert numpy.sqrt(rmsd_1 / 4) < 0.03
    assert numpy.sqrt(rmsd_2 / 4) < 0.03


def test_constants_h2o_ipoh():
    # Experimental data for validation is taken from http://www.ddbst.com/en/EED/VLE/VLE%202-Propanol%3BWater.php
    # Brunjes A.S.; Bogart M.J.P.: The Binary Systems Ethanol-n-Butanol, Acetone-Water
    # and Isopropanol-Water. Ind.Eng.Chem. 35 (1943) 255-260
    # NRTL Constants are taken from https://doi.org/10.1021/je960108n
    test_mixture = Mixtures.H2O_iPOH
    validation_compositions = [
        Composition(p=0.9796, type=CompositionType("molar")),
        Composition(p=0.7613, type=CompositionType("molar")),
        Composition(p=0.3029, type=CompositionType("molar")),
        Composition(p=0.06810, type=CompositionType("molar")),
    ]

    validation_temperature_list = [363.95, 354.26, 353.2, 354.36]

    validation_pressures = [
        (77.943036, 23.386964),
        (45.000653, 56.329347),
        (31.098177, 70.231823),
        (10.021537, 91.308463),
    ]
    tested_partial_pressures = [
        get_nrtl_partial_pressures(
            temperature=validation_temperature_list[i],
            mixture=test_mixture,
            composition=validation_compositions[i],
        )
        for i in range(4)
    ]

    rmsd_1 = 0
    rmsd_2 = 0

    for i in range(4):
        rmsd_1 = (
            tested_partial_pressures[i][0] - validation_pressures[i][0]
        ) ** 2 / validation_pressures[i][0] ** 2 + rmsd_1

        rmsd_2 = (
            tested_partial_pressures[i][1] - validation_pressures[i][1]
        ) ** 2 / validation_pressures[i][1] ** 2 + rmsd_2
        assert (
            abs(tested_partial_pressures[i][0] - validation_pressures[i][0])
            < validation_pressures[i][0] * 0.1
        )
        assert (
            abs(tested_partial_pressures[i][1] - validation_pressures[i][1])
            < validation_pressures[i][1] * 0.085
        )
    assert numpy.sqrt(rmsd_1 / 4) < 0.05
    assert numpy.sqrt(rmsd_2 / 4) < 0.05


def test_constants_etoh_etbe():
    # Experimental data for validation and NRTL constants are taken from
    # Isothermal vapor-liquid equilibria for binary and ternary systems containing ethyl tert-butyl ether,
    # ethanol, benzene, and toluene at 313.15 K
    # Oh, JH; Park, SJ
    # Journal of Industrial and Engineering Chemistry, 2005

    test_mixture = Mixtures.EtOH_ETBE
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
        for i in range(4)
    ]
    rmsd_1 = 0
    rmsd_2 = 0

    for i in range(4):
        rmsd_1 = (
            tested_partial_pressures[i][0] - validation_pressures[i][0]
        ) ** 2 / validation_pressures[i][0] ** 2 + rmsd_1

        rmsd_2 = (
            tested_partial_pressures[i][1] - validation_pressures[i][1]
        ) ** 2 / validation_pressures[i][1] ** 2 + rmsd_2
        assert (
            abs(tested_partial_pressures[i][0] - validation_pressures[i][0])
            < validation_pressures[i][0] * 0.05
        )
        assert (
            abs(tested_partial_pressures[i][1] - validation_pressures[i][1])
            < validation_pressures[i][1] * 0.05
        )
    assert numpy.sqrt(rmsd_1 / 4) < 0.03
    assert numpy.sqrt(rmsd_2 / 4) < 0.03


def test_constants_meoh_toluene():
    # Experimental data for validation is taken from https://doi.org/10.1016/0021-9614(88)90185-1
    # NRTL constants are taken from https://doi.org/10.1016/j.fluid.2019.112412

    test_mixture = Mixtures.MeOH_Toluene
    validation_compositions = [
        Composition(p=0.1830, type=CompositionType("molar")),
        Composition(p=0.4980, type=CompositionType("molar")),
        Composition(p=0.7640, type=CompositionType("molar")),
        Composition(p=0.960, type=CompositionType("molar")),
    ]

    validation_pressures = [
        (31.176927, 9.260073),
        (35.795292, 8.560708),
        (38.214708, 7.661292),
        (42.545305, 2.957695),
    ]

    tested_partial_pressures = [
        get_nrtl_partial_pressures(
            temperature=318,
            mixture=test_mixture,
            composition=validation_compositions[i],
        )
        for i in range(4)
    ]
    rmsd_1 = 0
    rmsd_2 = 0

    for i in range(4):
        rmsd_1 = (
            tested_partial_pressures[i][0] - validation_pressures[i][0]
        ) ** 2 / validation_pressures[i][0] ** 2 + rmsd_1

        rmsd_2 = (
            tested_partial_pressures[i][1] - validation_pressures[i][1]
        ) ** 2 / validation_pressures[i][1] ** 2 + rmsd_2
        assert (
            abs(tested_partial_pressures[i][0] - validation_pressures[i][0])
            < validation_pressures[i][0] * 0.05
        )
        assert (
            abs(tested_partial_pressures[i][1] - validation_pressures[i][1])
            < validation_pressures[i][1] * 0.05
        )
    assert numpy.sqrt(rmsd_1 / 4) < 0.03
    assert numpy.sqrt(rmsd_2 / 4) < 0.03


def test_constants_meoh_mtbe():
    # Experimental data for validation and NRTL Parameters are taken from
    # https://doi.org/10.1002/1521-4125(20020709)25:7<729::AID-CEAT729>3.0.CO;2-B

    test_mixture = Mixtures.MeOH_MTBE
    validation_compositions = [
        Composition(p=0.87425, type=CompositionType("molar")),
        Composition(p=0.59990, type=CompositionType("molar")),
        Composition(p=0.25010, type=CompositionType("molar")),
        Composition(p=0.0999, type=CompositionType("molar")),
    ]

    validation_pressures = [
        (61.746843, 31.823157),
        (41.161443, 52.408557),
        (25.82532, 67.74468),
        (14.484636, 79.085364),
    ]

    validation_temperatures = [328.6, 323.71, 322.71, 323.88]

    tested_partial_pressures = [
        get_nrtl_partial_pressures(
            temperature=validation_temperatures[i],
            mixture=test_mixture,
            composition=validation_compositions[i],
        )
        for i in range(4)
    ]
    rmsd_1 = 0
    rmsd_2 = 0

    for i in range(4):
        rmsd_1 = (
            tested_partial_pressures[i][0] - validation_pressures[i][0]
        ) ** 2 / validation_pressures[i][0] ** 2 + rmsd_1
        rmsd_2 = (
            tested_partial_pressures[i][1] - validation_pressures[i][1]
        ) ** 2 / validation_pressures[i][1] ** 2 + rmsd_2
        assert (
            abs(tested_partial_pressures[i][0] - validation_pressures[i][0])
            < validation_pressures[i][0] * 0.06
        )
        assert (
            abs(tested_partial_pressures[i][1] - validation_pressures[i][1])
            < validation_pressures[i][1] * 0.06
        )
    assert numpy.sqrt(rmsd_1 / 4) < 0.03
    assert numpy.sqrt(rmsd_2 / 4) < 0.03


def test_constants_meoh_dmc():
    # Experimental data for validation and NRTL Parameters are taken from https://doi.org/10.1016/j.fluid.2011.08.007

    test_mixture = Mixtures.MeOH_DMC
    validation_compositions = [
        Composition(p=0.096, type=CompositionType("molar")),
        Composition(p=0.318, type=CompositionType("molar")),
        Composition(p=0.814, type=CompositionType("molar")),
        Composition(p=0.905, type=CompositionType("molar")),
    ]

    validation_pressures = [
        (26.73066, 39.92934),
        (41.92914, 24.73086),
        (55.12782, 11.53218),
        (58.59414, 8.06586),
    ]

    validation_temperatures = [338.57, 330.24, 326.43, 326.52]

    tested_partial_pressures = [
        get_nrtl_partial_pressures(
            temperature=validation_temperatures[i],
            mixture=test_mixture,
            composition=validation_compositions[i],
        )
        for i in range(4)
    ]

    rmsd_1 = 0
    rmsd_2 = 0

    for i in range(4):
        rmsd_1 = (
            tested_partial_pressures[i][0] - validation_pressures[i][0]
        ) ** 2 / validation_pressures[i][0] ** 2 + rmsd_1
        rmsd_2 = (
            tested_partial_pressures[i][1] - validation_pressures[i][1]
        ) ** 2 / validation_pressures[i][1] ** 2 + rmsd_2
        assert (
            abs(tested_partial_pressures[i][0] - validation_pressures[i][0])
            < validation_pressures[i][0] * 0.06
        )
        assert (
            abs(tested_partial_pressures[i][1] - validation_pressures[i][1])
            < validation_pressures[i][1] * 0.06
        )
    assert numpy.sqrt(rmsd_1 / 4) < 0.04
    assert numpy.sqrt(rmsd_2 / 4) < 0.04
