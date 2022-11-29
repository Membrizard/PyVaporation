import typing
from copy import copy
from pathlib import Path

import attr
import joblib
import json
import numpy
from scipy import optimize

from ..diffusion_curve import DiffusionCurve, DiffusionCurveSet
from ..plotting import plot_graph, plot_surface


@attr.s(auto_attribs=True)
class Measurement:
    """
    Parameters:
    x (float): concentration
    t (float): temperature in K
    p (float): Permeance value
    """

    x: float
    t: float
    p: float


@attr.s(auto_attribs=True)
class Measurements:
    """
    List of Measurement objects
    """

    data: typing.List[Measurement]

    def __len__(self) -> int:
        return len(self.data)

    def __getitem__(self, item: int) -> Measurement:
        return self.data[item]

    def __add__(self, other):
        return Measurements(data=self.data + other.data)

    @classmethod
    def from_diffusion_curve_first(cls, curve: DiffusionCurve) -> "Measurements":
        """
        :param curve: Diffusion curve
        :return: Measurements for the first component
        """
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
        """
        :param curve: Diffusion curve
        :return: Measurements for the second component
        """
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
        """
        :param curves: Diffusion curves
        :return: Measurements for the first component from multiple diffusion curves
        """
        output = Measurements(data=[])
        for curve in curves:
            output += cls.from_diffusion_curve_first(curve)
        return output

    @classmethod
    def from_diffusion_curves_second(cls, curves: DiffusionCurveSet) -> "Measurements":
        """
        :param curves: Diffusion curves
        :return: Measurements for the second component from multiple diffusion curves
        """
        output = Measurements(data=[])
        for curve in curves:
            output += cls.from_diffusion_curve_second(curve)
        return output


@attr.s(auto_attribs=True)
class PervaporationFunction:
    """
    :param n: max power of temperature independent components in equation
    :param m: max power of temperature dependent components in equation
    :param alpha: multiplication coefficient
    :param a: list of coefficients for temperature independent components
    :param b: list of coefficients for temperature dependent components
    """

    n: int
    m: int

    alpha: float
    a: typing.List[float]
    b: typing.List[float]

    @classmethod
    def from_array(
        cls, array: typing.Union[typing.List[float], numpy.ndarray], n: int, m: int
    ) -> "PervaporationFunction":
        """

        :param array: list of coefficients
        :param n: max power of temperature independent components in equation
        :param m: max power of temperature dependent components in equation
        :return: PervaporationFunction
        """
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

    @classmethod
    def load(cls, path: typing.Union[str, Path]) -> "PervaporationFunction":
        """
        :param path: path to binary file with PervaporationFunction
        :return: PervaporationFunction
        """
        if type(path) is not Path:
            path = Path(path)
        return joblib.load(path)

    def save(self, path: typing.Union[str, Path]) -> None:
        """
        :param path: path where to save binary file with PervaporationFunction
        """
        joblib.dump(self, path)

    @classmethod
    def safe_load(cls, path: typing.Union[str, Path]) -> "PervaporationFunction":
        """
        :param path: Path to a json object
        :return: PervaporationFunction from a json file
        """
        with open(path, "r") as openfile:
            # Reading from json file
            json_object = json.load(openfile)
        return PervaporationFunction(
            n=json_object["n"],
            m=json_object["m"],
            alpha=json_object["alpha"],
            a=json_object["a"],
            b=json_object["b"],
        )

    def safe_save(self, path: typing.Union[str, Path]):
        """
        :param path: Path to a json object
        :return: Saves a PervaporationFunction to a json file
        """
        json_dict = {
            "n": self.n,
            "m": self.m,
            "alpha": self.alpha,
            "a": list(self.a),
            "b": list(self.b),
        }
        with open(path, "w") as outfile:
            json.dump(json_dict, outfile)

    def plot(
        self,
        experimental_data: typing.Optional[Measurements] = None,
        temperature: typing.Optional[float] = None,
        concentration: typing.Tuple[float, float] = None,
    ):
        points = {}
        x, t, p = 0, 0, 0

        if experimental_data is not None:
            x = [m.x for m in experimental_data]
            t = [m.t for m in experimental_data]
            p = [m.p for m in experimental_data]
            x_max = max(x)
            x_min = min(x)
            t_max = max(t)
            t_min = min(t)
            if temperature is not None and not (temperature in set(t)):
                raise ValueError(f"No experimental points at {temperature} K available")
            if temperature is not None:
                x_selected = []
                p_selected = []
                for i in range(len(t)):
                    if t[i] == temperature:
                        x_selected.append(x[i])
                        p_selected.append(p[i])
                x = x_selected
                p = p_selected

            points["Experiment"] = (x, p, False)
            if len(set(t)) == 1:
                temperature = set(t).pop()

        else:
            if concentration is None:
                x_max = 1
                x_min = 0
                t_max = 373.15
                t_min = 293.15
            else:
                if max(concentration) > 1 or min(concentration) < 0:
                    raise ValueError("Concentration must be in [0,1] range")
                x_max = max(concentration)
                x_min = min(concentration)
                t_max = 373.15
                t_min = 293.15

        x_v = numpy.linspace(x_min, x_max, num=50)

        if temperature is not None:
            p_fit = [self(x_v[i], temperature) for i in range(len(x_v))]
            points["Fit"] = (x_v, p_fit, True)
            plot_graph(
                x_label="First component fraction in feed",
                y_label="Permeance",
                points=points,
                title=f"Temperature {temperature} K",
            )
            return

        plot_surface(
            condition=experimental_data is not None,
            function=self,
            x=x,
            t=t,
            p=p,
            t_min=t_min,
            t_max=t_max,
            x_v=x_v,
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
    """
    Fit PervaporationFunction to given data with pre-defined parameters
    :param data: List of Measurements
    :param n: max power of temperature independent components in equation
    :param m: max power of temperature dependent components in equation
    :param include_zero: include zero-datapoint or not
    :param component_index: index of component to fit
    :return: PervaporationFunction
    """
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
    """
    Finds best fit of PervaporationFunction to given data
    :param data: List of Measurements
    :param include_zero: include zero-datapoint or not
    :param component_index: index of component to fit
    :param n: max power of temperature independent components in equation
    :param m: max power of temperature dependent components in equation
    :return: PervaporationFunction
    """
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
