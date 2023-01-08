from pyvaporation.mixtures.UNIQUAC_fitting import fit, VLEPoints

points = VLEPoints.from_csv(path="tests/VLE_data/H2O_EtOH.csv")

params = fit(data=points)

print(params)
