import typing
import attr
import numpy

from copy import copy
from ..utils import UNIQUACParameters
from ..components import Component
from .mixture import get_partial_pressures, Composition, Mixture
from scipy import optimize


@attr.s(auto_attribs=True)
class VLEPoint:
    composition: Composition
    pressures: typing.Tuple[float, float]
    temperature: float


@attr.s(auto_attribs=True)
class VLEPoints:
    """
        Class to represent experimental VLE points
    """
    data: typing.List[VLEPoint]

    @classmethod
    def from_file(
            cls, path: str,
    ) -> "VLEPoints":
        assert len(array) == 5
        return cls(
            alpha_12=array[0],
            alpha_21=array[1],
            beta_12=array[2],
            beta_21=array[3],
            z=int(array[4]),
        )

    def __len__(self):
        return len(self.data)

    def __getitem__(self, item: int) -> VLEPoint:
        return self.data[item]

    def __add__(self, other):
        return VLEPoints(data=self.data + other.data)


def fit(
    data: VLEPoints,
    components: typing.List[Component]
) -> UNIQUACParameters:
    """
    Get the UNIQUAC parameters to fit the VLE of a given mixture
    """
    _data = copy(data)
    result = optimize.minimize(
        lambda params: objective(data=_data,
                                 components=components,
                                 params=params),
        x0=numpy.array([0, 0, 0, 0, 10]),
        method="Powell",
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
