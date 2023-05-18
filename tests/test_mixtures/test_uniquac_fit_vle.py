from pyvaporation.activity_coefficient_models.uniquac_fitting import fit_vle, VLEPoints
from pyvaporation.mixtures import get_partial_pressures, BinaryMixture


def test_fit_uniquac_vle():

    points = VLEPoints.from_csv(path=f"tests/VLE_data/binary/H2O_EtOH.csv")

    params = fit_vle(data=points, method="COBYLA")

    test_mixture = BinaryMixture(
        name="",
        first_component=points.components[0],
        second_component=points.components[1],
        uniquac_params=params,
    ).to_mixture()

    errors_h2o = []
    errors_etoh = []

    for point in points:
        calc_pressures = get_partial_pressures(
            temperature=point.temperature,
            mixture=test_mixture,
            composition=point.composition,
            calculation_type="UNIQUAC",
        )

        errors_h2o.append(
            abs((point.pressures[0] - calc_pressures[0]) / point.pressures[0])
        )
        errors_etoh.append(
            abs((point.pressures[1] - calc_pressures[1]) / point.pressures[1])
        )

    assert max(errors_h2o) < 0.1
    assert max(errors_etoh) < 0.1


# TODO Add tests for multicomponent mixtures fitting, reading and other
