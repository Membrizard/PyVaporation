from components import Components


def test_constants_h2o():

    # Antoine Pressure Constants Calculated from https://dx.doi.org/10.1115/1.3687121

    assert abs(Components.H2O.get_vapor_pressure(313) - 7.31934) < 1e-4
    assert abs(Components.H2O.get_vapor_pressure(323) - 12.24821) < 1e-4
    assert abs(Components.H2O.get_vapor_pressure(333) - 19.78961) < 1e-4

    # Specific Heat Constants Calculated from http://dx.doi.org/10.1021/je900208n

    assert abs(Components.H2O.get_specific_heat(313) - 33.725516) < 1e-3
    assert abs(Components.H2O.get_specific_heat(323) - 33.800943) < 1e-3
    assert abs(Components.H2O.get_specific_heat(333) - 33.877785) < 1e-3


def test_constants_meoh():

    # Antoine Pressure Constants Calculated from https://dx.doi.org/10.1016/0021-9614(75)90267-0

    assert abs(Components.MeOH.get_vapor_pressure(313) - 35.1986) < 1e-4
    assert abs(Components.MeOH.get_vapor_pressure(323) - 55.2130) < 1e-4
    assert abs(Components.MeOH.get_vapor_pressure(333) - 84.0491) < 1e-4

    # Specific Heat Constants Calculated from https://doi.org/10.1007/BF00870444

    assert abs(Components.MeOH.get_specific_heat(313) - 84.24440) < 1e-3
    assert abs(Components.MeOH.get_specific_heat(323) - 86.61974) < 1e-3
    assert abs(Components.MeOH.get_specific_heat(333) - 89.21497) < 1e-3


def test_constants_etoh():

    # Antoine Pressure Constants Calculated from https://dx.doi.org/10.1016/0021-9614(75)90267-0

    assert abs(Components.EtOH.get_vapor_pressure(313) - 17.7708) < 1e-4
    assert abs(Components.EtOH.get_vapor_pressure(323) - 29.2779) < 1e-4
    assert abs(Components.EtOH.get_vapor_pressure(333) - 46.5843) < 1e-4

    # Specific Heat Constants Calculated from Mazur, V.J., On the specific heat of ethyl alcohol,
    # Acta Phys. Polon., 1940, 8, 6-11.

    assert abs(Components.EtOH.get_specific_heat(313) - 122.1500349) < 1e-3
    assert abs(Components.EtOH.get_specific_heat(323) - 127.4315734) < 1e-3
    assert abs(Components.EtOH.get_specific_heat(333) - 133.0909154) < 1e-3


def test_constants_ipoh():

    # Antoine Pressure Constants Calculated from https://doi.org/10.1016/j.fluid.2015.09.052

    assert abs(Components.iPOH.get_vapor_pressure(313) - 13.7388) < 1e-4
    assert abs(Components.iPOH.get_vapor_pressure(323) - 23.4142) < 1e-4
    assert abs(Components.iPOH.get_vapor_pressure(333) - 38.2827) < 1e-4

    # Specific Heat Constants Calculated from https://doi.org/10.1021/je60083a031

    assert abs(Components.iPOH.get_specific_heat(313) - 167.52137) < 1e-3
    assert abs(Components.iPOH.get_specific_heat(323) - 177.35520) < 1e-3
    assert abs(Components.iPOH.get_specific_heat(333) - 188.25157) < 1e-3


def test_constants_dme():

    # Antoine Pressure Constants Calculated from https://dx.doi.org/10.1021/ie50448a022

    assert abs(Components.DME.get_vapor_pressure(313) - 18.4210) < 1e-4
    assert abs(Components.DME.get_vapor_pressure(323) - 26.6299) < 1e-4
    assert abs(Components.DME.get_vapor_pressure(333) - 37.5490) < 1e-4

    # Specific Heat Constants Calculated from http://doi.org/10.1021/acs.jced.1c00229

    assert abs(Components.DME.get_specific_heat(313) - 193.95434) < 1e-3
    assert abs(Components.DME.get_specific_heat(323) - 195.75853) < 1e-3
    assert abs(Components.DME.get_specific_heat(333) - 197.73868) < 1e-3


def test_constants_dmc():

    # Antoine Pressure Constants Calculated from https://doi.org/10.1016/j.fluid.2011.08.007

    assert abs(Components.DMC.get_vapor_pressure(313) - 14.6761) < 1e-4
    assert abs(Components.DMC.get_vapor_pressure(323) - 22.9947) < 1e-4
    assert abs(Components.DMC.get_vapor_pressure(333) - 34.7419) < 1e-4

    # Specific Heat Constants Calculated from http://doi.org/10.1021/acs.jced.7b00295

    assert abs(Components.DMC.get_specific_heat(313) - 167.09697) < 1e-3
    assert abs(Components.DMC.get_specific_heat(323) - 168.69313) < 1e-3
    assert abs(Components.DMC.get_specific_heat(333) - 170.46484) < 1e-3


def test_constants_mtbe():

    # Antoine Pressure Constants Calculated from https://dx.doi.org/10.1135/cccc19691317

    assert abs(Components.MTBE.get_vapor_pressure(313) - 60.1905) < 1e-4
    assert abs(Components.MTBE.get_vapor_pressure(323) - 85.8647) < 1e-4
    assert abs(Components.MTBE.get_vapor_pressure(333) - 119.4935) < 1e-4

    # Specific Heat Constants Calculated from https://doi.org/10.1016/0021-9614(75)90194-9

    assert abs(Components.MTBE.get_specific_heat(313) - 192.4237) < 1e-3
    assert abs(Components.MTBE.get_specific_heat(323) - 195.8651) < 1e-3
    assert abs(Components.MTBE.get_specific_heat(333) - 199.3988) < 1e-3


def test_constants_etbe():

    # Antoine Pressure Constants Calculated from https://dx.doi.org/10.1135/cccc19691317

    assert abs(Components.ETBE.get_vapor_pressure(313) - 31.2979) < 1e-4
    assert abs(Components.ETBE.get_vapor_pressure(323) - 45.9826) < 1e-4
    assert abs(Components.ETBE.get_vapor_pressure(333) - 65.7505) < 1e-4

    # Specific Heat Constants Calculated from https://doi.org/10.1021/je900208n

    assert abs(Components.ETBE.get_specific_heat(313) - 226.5025) < 1e-3
    assert abs(Components.ETBE.get_specific_heat(323) - 230.8717) < 1e-3
    assert abs(Components.ETBE.get_specific_heat(333) - 235.3427) < 1e-3


def test_constants_cyclohexane():

    # Antoine Pressure Constants Calculated from https://dx.doi.org/10.1039/tf9686400637

    assert abs(Components.CycloHexane.get_vapor_pressure(313) - 24.6260) < 1e-4
    assert abs(Components.CycloHexane.get_vapor_pressure(323) - 36.2458) < 1e-4
    assert abs(Components.CycloHexane.get_vapor_pressure(333) - 51.9116) < 1e-4

    # Specific Heat Constants Calculated from Safir, L.I., Experimental determination
    # of the isobaric heat capacity of cyclohexane at atmospheric pressure,
    # Izv. Vyssh. Uchebn. Zaved. Neft. Gaz 21, 1978, (12), 81-82.

    assert abs(Components.CycloHexane.get_specific_heat(313) - 162.4372) < 1e-3
    assert abs(Components.CycloHexane.get_specific_heat(323) - 166.2413) < 1e-3
    assert abs(Components.CycloHexane.get_specific_heat(333) - 170.0453) < 1e-3


def test_constants_benzene():

    # Antoine Pressure Constants Calculated from https://doi.org/10.1002/9781118477304.app2

    assert abs(Components.Benzene.get_vapor_pressure(313) - 24.3752) < 1e-4
    assert abs(Components.Benzene.get_vapor_pressure(323) - 36.1835) < 1e-4
    assert abs(Components.Benzene.get_vapor_pressure(333) - 52.2136) < 1e-4

    # Specific Heat Constants Calculated from https://doi.org/10.1007/BF00503954

    assert abs(Components.Benzene.get_specific_heat(313) - 139.9297) < 1e-3
    assert abs(Components.Benzene.get_specific_heat(323) - 142.3563) < 1e-3
    assert abs(Components.Benzene.get_specific_heat(333) - 144.8910) < 1e-3


def test_constants_toluene():

    # Antoine Pressure Constants Calculated from https://doi.org/10.1002/9781118477304.app2

    assert abs(Components.Toluene.get_vapor_pressure(313) - 7.8949) < 1e-4
    assert abs(Components.Toluene.get_vapor_pressure(323) - 12.2889) < 1e-4
    assert abs(Components.Toluene.get_vapor_pressure(333) - 18.5332) < 1e-4

    # Specific Heat Constants Calculated from https://doi.org/10.1021/j100573a011

    assert abs(Components.Toluene.get_specific_heat(313) - 161.09761) < 1e-3
    assert abs(Components.Toluene.get_specific_heat(323) - 163.93334) < 1e-3
    assert abs(Components.Toluene.get_specific_heat(333) - 166.80588) < 1e-3


def test_constants_acetic_acid():

    # Antoine Pressure Constants Calculated from https://dx.doi.org/10.1021/je60004a009

    assert abs(Components.AceticAcid.get_vapor_pressure(313) - 4.68410) < 1e-4
    assert abs(Components.AceticAcid.get_vapor_pressure(323) - 7.63603) < 1e-4
    assert abs(Components.AceticAcid.get_vapor_pressure(333) - 12.0401) < 1e-4

    # Specific Heat Constants Calculated from
    # Parks G.S.; Kelley K.K.: Thermal Data on Organic Compounds. II. The Heat Capacities of Five Organic Compounds.
    # The Entropies and Free Energies of Some Homologous Series of Aliphatic Compounds.
    # J.Am.Chem.Soc. 47 (1925) 2089-2097
    # and
    # von Reis M.A.: Die spezifische Wärme flüssiger organischer Verbindungen und
    # ihre Beziehung zu deren Molekulargewicht. Ann.Physik 249 (1881) 447-465

    assert abs(Components.AceticAcid.get_specific_heat(292.6) - 122.61) < 6e-1
    assert abs(Components.AceticAcid.get_specific_heat(334.45) - 128.34) < 6e-1
    assert abs(Components.AceticAcid.get_specific_heat(358.55) - 130.80) < 6e-1

