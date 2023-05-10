from pyvaporation.mixtures import Mixtures, get_partial_pressures, Composition
from pyvaporation.components import Components
from pyvaporation.mixtures import VLEPoints, fit_vle, Mixture
import numpy

ternary_points = VLEPoints.from_csv(path="./tests/VLE_data/ternary/H2O_MeOH_EtOH.csv")
binary_points = VLEPoints.from_csv(path="./tests/VLE_data/binary/H2O_EtOH.csv")
ternary_coefs = fit_vle(data=ternary_points)
binary_coefs = fit_vle(data=binary_points)
print(f"Binary: {binary_coefs} \nTernary: {ternary_coefs} ")

test_mixture_1 = Mixture(
    name="Ternary",
    components=[Components.MeOH, Components.EtOH, Components.H2O],
    uniquac_params=ternary_coefs,
)


error_1 = 0
error_2 = 0
error_3 = 0
validation_average_1 = 0
validation_average_2 = 0
validation_average_3 = 0

sum_error_1 = []
sum_error_2 = []
sum_error_3 = []


for point in ternary_points:
    calc_pressures = get_partial_pressures(
            temperature=point.temperature,
            mixture=test_mixture_1,
            composition=point.composition,
            calculation_type="UNIQUAC",
        )
    print(calc_pressures)

    validation_average_1 += point.pressures[0]
    validation_average_2 += point.pressures[1]
    validation_average_3 += point.pressures[2]

    error_1 += (calc_pressures[0] - point.pressures[0]) ** 2
    error_2 += (calc_pressures[1] - point.pressures[1]) ** 2
    error_3 += (calc_pressures[2] - point.pressures[2]) ** 2

    sum_error_1.append(numpy.sqrt(error_1 * len(ternary_points)) / validation_average_1)
    sum_error_2.append(numpy.sqrt(error_2 * len(ternary_points)) / validation_average_2)
    sum_error_3.append(numpy.sqrt(error_3 * len(ternary_points)) / validation_average_3)

    print(f"Ternary:\nError 1: {max(sum_error_1)}")
    print(f"Error 2: {max(sum_error_2)}")
    print(f"Error 3: {max(sum_error_3)}")

test_mixture_2 = Mixture(
    name="Binary",
    components=[Components.H2O, Components.EtOH],
    uniquac_params=binary_coefs,
)


error_1 = 0
error_2 = 0
validation_average_1 = 0
validation_average_2 = 0
sum_error_1 = []
sum_error_2 = []

for point in binary_points:
    calc_pressures = get_partial_pressures(
            temperature=point.temperature,
            mixture=test_mixture_2,
            composition=point.composition,
            calculation_type="UNIQUAC",
        )
    print(calc_pressures)

    validation_average_1 += point.pressures[0]
    validation_average_2 += point.pressures[1]

    error_1 += (calc_pressures[0] - point.pressures[0]) ** 2
    error_2 += (calc_pressures[1] - point.pressures[1]) ** 2

    sum_error_1.append(numpy.sqrt(error_1 * len(binary_points)) / validation_average_1)
    sum_error_2.append(numpy.sqrt(error_2 * len(binary_points)) / validation_average_2)

    print(f"Binary:\nError 1: {max(sum_error_1)}")
    print(f"Error 2: {max(sum_error_2)}")
