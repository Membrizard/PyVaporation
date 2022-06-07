import typing
import numpy
import attr

from scipy import optimize


@attr.s(auto_attribs=True)
class Measurement:
    x: float
    t: float
    p: float


@attr.s(auto_attribs=True)
class FitData:
    measurements: typing.List[Measurement]

    def __len__(self) -> int:
        return len(self.measurements)

    def __getitem__(self, item: int) -> Measurement:
        return self.measurements[item]


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
        assert len(array) == 3 + n + m
        return cls(
            n=n,
            m=m,
            alpha=array[0],
            a=array[1: n + 2],
            b=array[n + 2:],
        )

    def __call__(self, x: float, t: float) -> float:
        return self.alpha * (
            numpy.exp(
                sum(self.a[i] * x**i for i in range(len(self.a)))
                - sum(self.b[i] * x**i for i in range(len(self.b))) / t
            )
        )


def get_initial_guess(n: int, m: int) -> typing.List[float]:
    return [1] * (3 + n + m)


def objective(data: FitData, params: typing.List[float], n: int, m: int) -> float:
    error = 0
    foo = PervaporationFunction.from_array(array=params, n=n, m=m)
    for d in data:
        error += (foo(d.x, d.t) - d.p) ** 2
    return numpy.sqrt(error / len(data))


def fit(data: FitData, n: int, m: int) -> PervaporationFunction:
    result = optimize.minimize(
        lambda params: objective(data, params, n=n, m=m),
        x0=numpy.array([1] * (3 + n + m)),
        method="Powell",
    )
    return PervaporationFunction.from_array(array=result.x, n=n, m=m)
