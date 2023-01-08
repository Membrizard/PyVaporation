import typing
import attr
import numpy
import pandas
from pathlib import Path

from copy import copy
from ..utils import UNIQUACParameters
from ..components import Component, Components
from .mixture import get_partial_pressures, Composition, Mixture
from scipy import optimize

VLE_COLUMNS = [
    'first_component',
    'second_component',
    'composition',
    'composition_type',
    'first_component_pressure',
    'second_component_pressure',
    'temperature',
    'reference']


@attr.s(auto_attribs=True)
class VLEPoint:
    """
            Class to represent single experimental VLE point
    """
    composition: Composition
    pressures: typing.Tuple[float, float]
    temperature: float

    @classmethod
    def from_dict(
            cls, d: typing.Mapping[str, typing.Union[str, float]]
    ) -> "VLEPoint":
        """
        Reads VLEPoint objects from dictionary
        :param d: Mapping of VLEPoint parameters
        :return: VLEPoint with parameters specified in mapping
        """

        composition = Composition(p=d["composition"], type=d["composition_type"])
        return cls(
            composition=composition,
            pressures=(d["first_component_pressure"], d["second_component_pressure"]),
            temperature=d["temperature"],
        )


@attr.s(auto_attribs=True)
class VLEPoints:
    """
        Class to represent experimental VLE points
    """
    components: typing.List[Component]
    data: typing.List[VLEPoint]

    @classmethod
    def from_csv(cls, path: typing.Union[str, Path]) -> "VLEPoints":
        """
        Generates VLEPOINTS object from the .csv file
        :param path: Path to the file
        :return: VLEPoints object with data described in a .csv file
        """
        frame = pandas.read_csv(path)

        if list(frame.columns) != VLE_COLUMNS:
            raise ValueError("Incorrect columns: %s" % list(frame.columns))

        points = []

        d = frame.iloc[0].to_dict()
        components = [getattr(Components, d["first_component"]),
                      getattr(Components, d["second_component"])]

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
            raise ValueError("Both sets should have same components and indexes")

        return VLEPoints(components=self.components, data=self.data + other.data)


def fit(
    data: VLEPoints,
) -> UNIQUACParameters:
    """
    Get the UNIQUAC parameters to fit the VLE of a given mixture
    """
    _data = copy(data)
    result = optimize.minimize(
        lambda params: objective(data=_data,
                                 components=data.components,
                                 params=params),
        x0=numpy.array([0, 0, 0, 0, 10]),
        method="Nelder-Mead",
    )
    return UNIQUACParameters.from_array(result.x)


def objective(data: VLEPoints, components: typing.List[Component], params: typing.List[float]) -> float:
    error = 0
    mixture = Mixture(
                      name="",
                      first_component=components[0],
                      second_component=components[1],
                      uniquac_params=UNIQUACParameters.from_array(params),
    )
    for point in data:
        error += (get_partial_pressures(temperature=point.temperature,
                                        composition=point.composition,
                                        mixture=mixture,
                                        calculation_type="UNIQUAC",
                                        )[0] - point.pressures[0]) ** 2 \
                 + (get_partial_pressures(temperature=point.temperature,
                                          composition=point.composition,
                                          mixture=mixture,
                                          calculation_type="UNIQUAC",
                                          )[1] - point.pressures[1]) ** 2
    return numpy.sqrt(error / len(data))
