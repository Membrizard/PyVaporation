from pyvaporation.mixtures.uniquac_fitting import fit_vle, VLEPoints


def test_fit_uniquac_vle():
    # TODO: fix test in pipeline
    points = VLEPoints.from_csv(path=f"tests/VLE_data/binary/H2O_EtOH.csv")

    params = fit_vle(data=points, method=None)

    assert abs(params.alpha_12-21.127561704493143) < 1e-5
    assert abs(params.alpha_21 - 100.10268878024358) < 1e-5
    assert abs(params.beta_12 - -0.9175664931087569) < 1e-5
    assert abs(params.beta_21 - 2.4619377106475753) < 1e-5
    assert params.z == 13



