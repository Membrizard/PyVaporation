import typing
from copy import copy

import attr
import numpy
from scipy import optimize

from diffusion_curve import DiffusionCurve, DiffusionCurveSet


@attr.s(auto_attribs=True)
class Measurement:
    x: float
    t: float
    p: float


@attr.s(auto_attribs=True)
class Measurements:
    data: typing.List[Measurement]

    def __len__(self) -> int:
        return len(self.data)

    def __getitem__(self, item: int) -> Measurement:
        return self.data[item]

    def __add__(self, other):
        return Measurements(data=self.data + other.data)

    @classmethod
    def from_diffusion_curve_first(cls, curve: DiffusionCurve) -> "Measurements":
        return Measurements(
            data=[
                Measurement(
                    x=curve.feed_compositions[i].first,
                    t=curve.feed_temperature,
                    p=curve.permeances[i][0].value,
                )
                for i in range(len(curve))
            ]
        )

    def append(self, measurement: Measurement) -> None:
        self.data.append(measurement)

    @classmethod
    def from_diffusion_curve_second(cls, curve: DiffusionCurve) -> "Measurements":
        return Measurements(
            data=[
                Measurement(
                    x=curve.feed_compositions[i].first,
                    t=curve.feed_temperature,
                    p=curve.permeances[i][1].value,
                )
                for i in range(len(curve))
            ]
        )

    @classmethod
    def from_diffusion_curves_first(cls, curves: DiffusionCurveSet) -> "Measurements":
        output = Measurements(data=[])
        for curve in curves:
            output += cls.from_diffusion_curve_first(curve)
        return output

    @classmethod
    def from_diffusion_curves_second(cls, curves: DiffusionCurveSet) -> "Measurements":
        output = Measurements(data=[])
        for curve in curves:
            output += cls.from_diffusion_curve_second(curve)
        return output


@attr.s(auto_attribs=True)
class PervaporationFunction:
    n: int
    m: int

    alpha: float
    a: typing.List[float]
    b: typing.List[float]

    @classmethod
    def from_array(
        cls, array: typing.Union[typing.List[float], numpy.ndarray], n: int, m: int
    ) -> "PervaporationFunction":
        assert len(array) == 2 + n + m
        return cls(
            n=n,
            m=m,
            alpha=array[0],
            a=array[1 : n + 1],
            b=array[n + 1 :],
        )

    def __call__(self, x: float, t: float) -> float:
        return self.alpha * (
            numpy.exp(
                sum(self.a[i] * numpy.power(x, i + 1) for i in range(len(self.a)))
                - sum(self.b[i] * numpy.power(x, i) for i in range(len(self.b))) / t
            )
        )

    def __mul__(self, constant: float) -> "PervaporationFunction":
        return PervaporationFunction(
            n=self.n,
            m=self.m,
            alpha=self.alpha * constant,
            a=self.a,
            b=self.b,
        )


def _suggest_n_m(
    data: Measurements, n: typing.Optional[int] = None, m: typing.Optional[int] = None
) -> typing.Tuple[int, int]:
    if n is None:
        suggested_n = int(numpy.sqrt(len(data)))
    else:
        suggested_n = n
    if m is None:
        suggested_m = int(numpy.sqrt(len(data)))
    else:
        suggested_m = m
    return suggested_n, suggested_m


def get_initial_guess(n: int, m: int) -> typing.List[float]:
    return [1] * (2 + n + m)


def objective(data: Measurements, params: typing.List[float], n: int, m: int) -> float:
    error = 0
    foo = PervaporationFunction.from_array(array=params, n=n, m=m)
    for d in data:
        error += (foo(d.x, d.t) - d.p) ** 2
    return numpy.sqrt(error / len(data))


def fit(
    data: Measurements,
    n: typing.Optional[int] = None,
    m: typing.Optional[int] = None,
    include_zero: bool = False,
    component_index: int = 0,
) -> PervaporationFunction:

    _n, _m = _suggest_n_m(data, n, m)
    if component_index not in {0, 1}:
        raise ValueError("Index should be either 0 or 1")

    _data = copy(data)

    if include_zero:
        unique_temperatures = set([m.t for m in data])
        for t in unique_temperatures:
            _data.append(
                Measurement(
                    x=component_index,
                    t=t,
                    p=0,
                )
            )

    result = optimize.minimize(
        lambda params: objective(_data, params, n=_n, m=_m),
        x0=numpy.array([0] * (2 + _n + _m)),
        method="Powell",
    )
    return PervaporationFunction.from_array(array=result.x, n=_n, m=_m)


def find_best_fit(
    data: Measurements,
    include_zero: bool = False,
    component_index: int = 0,
    n: typing.Optional[int] = None,
    m: typing.Optional[int] = None,
) -> PervaporationFunction:

    max_power_n = min(5, round(numpy.power(len(data), 0.5)))
    max_power_m = min(5, max(1, len(set([measurement.t for measurement in data])) - 1))

    if n is None:
        n_tries = list(range(max_power_n))
    else:
        n_tries = list(range(n + 1))
        if n >= len(data):
            print(
                "n is more or equal to the number of available points, you may be over-fitting"
            )
    if m is None:
        m_tries = list(range(max_power_m))
    else:
        m_tries = list(range(m + 1))
        if m >= len(set([measurement.t for measurement in data])):
            print(
                "m is more or equal to the number of available temperature points, you may be over-fitting"
            )

    best_curve = None
    best_loss = numpy.inf

    for n in n_tries:
        for m in m_tries:
            curve = fit(
                data,
                n=n,
                m=m,
                include_zero=include_zero,
                component_index=component_index,
            )
            loss = sum(
                [
                    (curve(measurement.x, measurement.t) - measurement.p) ** 2
                    for measurement in data
                ]
            )

            if loss < best_loss:
                best_curve = curve
                best_loss = loss

    return best_curve
