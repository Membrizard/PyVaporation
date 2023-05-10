from pyvaporation.mixtures import VLEPoints, fit_vle

points = VLEPoints.from_csv(path="tests/VLE_data/binary/MeOH_DMC.csv")
params = fit_vle(data=points)
print(params)
