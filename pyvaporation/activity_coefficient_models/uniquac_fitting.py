import math
import typing
import attr
import numpy
import pandas
from pathlib import Path


from ..utils import UNIQUACParameters
from ..components import Component, Components
from ..mixtures import get_partial_pressures, Composition, Mixture
from scipy import optimize

VLE_COLUMNS = [
    "first_component",
    "second_component",
    "composition",
    "composition_type",
    "first_component_pressure",
    "second_component_pressure",
    "temperature",
    "reference",
]

# VLE_COLUMNS = [
#     "components",
#     "composition",
#     "composition_type",
#     "first_component_pressure",
#     "second_component_pressure",
#     "temperature",
#     "reference",
# ]

FITTING_ALGS = [
    "Nelder-Mead",
    "Powell",
    "CG",
    "BFGS",
    "L-BFGS-B",
    "TNC",
    "COBYLA",
    "SLSQP",
    "trust-constr",
]


@attr.s(auto_attribs=True)
class VLEPoint:
    """
    Class to represent single experimental VLE point
    """

    composition: Composition
    pressures: typing.Tuple[float]
    temperature: float
    reference: typing.Optional[str]

    @classmethod
    def from_dict(cls, d: typing.Mapping[str, typing.Union[str, float]]) -> "VLEPoint":
        """
        Reads VLEPoint objects from dictionary
        :param d: Mapping of VLEPoint parameters
        :return: VLEPoint with parameters specified in mapping
        """

        if isinstance(d["composition"], str):
            p = [float(x) for x in d["composition"].split(sep=" ")]
        elif isinstance(d["composition"], float):
            p = d["composition"]
        else:
            raise ValueError("Unrecognisable type of Composition input.")

        pressures = tuple(float(d[x]) for x in list(d.keys()) if x.endswith("_component_pressure"))
        composition = Composition(p=p, type=d["composition_type"])

        if len(composition) != len(pressures):
            raise ValueError("The number of values for composition and components' pressures must correspond")

        return cls(
            composition=composition,
            pressures=pressures,
            temperature=d["temperature"],
            reference=d["reference"]
        )


@attr.s(auto_attribs=True)
class VLEPoints:
    """
    Class to represent experimental VLE points
    """

    components: typing.List[Component]
    data: typing.List[VLEPoint]

    def __attrs_post_init__(self):
        for value in self.data:
            if len(value.composition) != len(self.components):
                raise ValueError(f"Composition does not correspond to the list of components for {value}")
            if len(value.pressures) != len(self.components):
                raise ValueError(f"List of pressures does not correspond to the list of components for {value}")

    @classmethod
    def from_csv(cls, path: typing.Union[str, Path]) -> "VLEPoints":
        """
        Generates VLEPOINTS object from the .csv file
        :param path: Path to the file
        :return: VLEPoints object with data described in a .csv file
        """
        frame = pandas.read_csv(path)

        component_columns = [x for x in list(frame.columns) if x.endswith("_component")]

        if len(component_columns) < 2:
            raise ValueError("At least two components must be indicated.")

        points = []

        d = frame.iloc[0].to_dict()

        components = [getattr(Components, d[i]) for i in component_columns]

        for _, row in frame.iterrows():
            d = row.to_dict()
            points.append(VLEPoint.from_dict(d))

        return VLEPoints(components=components, data=points)

    def __len__(self):
        return len(self.data)

    def __getitem__(self, item: int) -> VLEPoint:
        return self.data[item]

    def __add__(self, other):
        if self.components != other.components:
            raise ValueError("Both sets should have same components in the same order")

        return VLEPoints(components=self.components, data=self.data + other.data)


def fit_vle(
    data: VLEPoints,
    method: typing.Optional[str] = None,
) -> UNIQUACParameters:
    """
    Get the UNIQUAC parameters to fit the VLE of a given mixture
    :param data Experimental data represented as a VLEPoints object
    :param method Optimization method string, if left None the most accurate method is chosen
    :return fitted UNIQUACParameters object
    """
    if method is None:
        algs = FITTING_ALGS
    else:
        algs = [method]

    best_fit = []
    error = 1000
    initial_guess = numpy.append(numpy.zeros(math.comb(len(data.components), 2)*4), 10)

    for alg in algs:
        result = optimize.minimize(
            lambda params: objective(data=data, params=params),
            x0=numpy.array(initial_guess),
            method=alg,
        )

        current_error = objective(data=data, params=result.x)
        if current_error < error:
            best_fit = result.x
            error = current_error

    return UNIQUACParameters.from_array(array=best_fit, components=[x.name for x in data.components])


def objective(data: VLEPoints, params: typing.List[float]) -> float:
    """
    Objective function for minimization during the UNIQUAC parameters fit
    :param data: data Experimental data represented as a VLEPoints object
    :param params: UNIQUAC parameters represented as a list
    :return: accumulative squared error as float
    """
    error = 0

    mixture = Mixture(
        name="fitting_mixture",
        components=data.components,
        uniquac_params=UNIQUACParameters.from_array(array=params, components=[i.name for i in data.components])
    )

    for point in data:
        calc = get_partial_pressures(
            temperature=point.temperature,
            composition=point.composition,
            mixture=mixture,
            calculation_type="UNIQUAC",
        )
        for i in range(len(point.pressures)):
            error += (calc[i]-point.pressures[i])**2

    return numpy.sqrt(error / len(data))
