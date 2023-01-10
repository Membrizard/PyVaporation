from pyvaporation.mixtures.uniquac_fitting import fit_vle, VLEPoints
from pyvaporation.mixtures import get_partial_pressures, Mixture
from pyvaporation.components import Components

files = ["H2O_EtOH.csv",
         "H2O_MeOH.csv",
         "H2O_iPOH.csv",
         "H2O_AceticAcid.csv",
         "EtOH_ETBE.csv",
         "MeOH_Toluene.csv",
         "MeOH_MTBE.csv",
         "MeOH_DMC.csv",
         ]




for file in files:
    points = VLEPoints.from_csv(path=f"tests/VLE_data/binary/{file}")

    params = fit_vle(data=points, method=None)

    print(f"parameters for {file} {params}")

    # Checking the accuracy

    test_mixture = Mixture(
        name="",
        first_component=points.components[0],
        second_component=points.components[1],
        uniquac_params=params
    )

    errors_h2o = []
    errors_etoh = []

    fc_calc = []
    sc_calc = []

    for point in points:
        calc_pressures = get_partial_pressures(temperature=point.temperature,
                                               mixture=test_mixture,
                                               composition=point.composition,
                                               calculation_type="UNIQUAC")

        #print(f"calculated:{calc_pressures} experimental:{point.pressures} error fc {(point.pressures[0]-calc_pressures[0])/point.pressures[0]} error sc {(point.pressures[1]-calc_pressures[1])/point.pressures[1]}")
        errors_h2o.append(abs((point.pressures[0]-calc_pressures[0])/point.pressures[0]))
        errors_etoh.append(abs((point.pressures[1]-calc_pressures[1])/point.pressures[1]))

        fc_calc.append(calc_pressures[0])
        sc_calc.append(calc_pressures[1])

    print(max(errors_h2o))
    print(max(errors_etoh))
    print(fc_calc)
    print(sc_calc)


