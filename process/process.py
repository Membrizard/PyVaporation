import typing
from pathlib import Path

import attr
import numpy
import pandas
import joblib

from random import randrange
from datetime import datetime
from conditions import Conditions
from mixtures import Composition, Mixture, get_nrtl_partial_pressures
from optimizer import PervaporationFunction
from permeance import Permeance
from plotting import plot_graph

PROCESS_MODEL_COLUMNS = [
    "membrane_name",
    "mixture",
    "time",
    "feed_temperature",
    "permeate_temperature",
    "permeate_pressure",
    "composition",
    "composition_type",
    "permeate_composition",
    "permeate_composition_type" "partial_flux_1",
    "partial_flux_2",
    "permeance_1",
    "permeance_2",
    "units",
    "feed_evaporation_heat",
    "permeate_condensation_heat",
    "comment",
]


@attr.s(auto_attribs=True)
class ProcessModel:
    mixture: Mixture
    membrane_name: str
    feed_temperature: typing.List[float]
    feed_compositions: typing.List[Composition]
    permeate_composition: typing.List[Composition]
    permeate_temperature: typing.List[float]
    permeate_pressure: typing.List[float]
    feed_mass: typing.List[float]
    partial_fluxes: typing.List[typing.Tuple[float, float]]
    permeances: typing.List[typing.Tuple[Permeance, Permeance]]
    time: typing.List[float]
    feed_evaporation_heat: typing.List[float]
    permeate_condensation_heat: typing.List[typing.Optional[float]]
    initial_conditions: typing.Optional[Conditions] = None
    permeance_fits: typing.Optional[
        typing.Tuple[PervaporationFunction, PervaporationFunction]
    ] = None
    comments: typing.Optional[str] = None
    results_path: typing.Optional[Path] = None

    def __attrs_post_init__(self):

        for i in range(len(self.time)):
            if (
                self.permeate_pressure[i] is None
                and self.permeate_temperature[i] is not None
            ):
                self.permeate_pressure[i] = sum(
                    get_nrtl_partial_pressures(
                        self.permeate_temperature[i],
                        self.mixture,
                        self.permeate_composition[i],
                    )
                )

    def plot(self, y: typing.List, y_label: str = "", curve: bool = 1):
        """
        Draws basic plot of a specified parameter versus process time
        :param y: parameter
        :param y_label: Name of the axis
        :param curve: if True - connects points with a line, if False - draws raw points
        :return: plots a graph
        """

        x = self.time
        if type(y[0]) == tuple:
            first = []
            second = []
            if isinstance(y[0][0], Permeance):
                for m in y:
                    first.append(m[0].value)
                    second.append(m[1].value)
            elif type(y[0][0]) == float or type(y[0][0]) == numpy.float64:
                for m in y:
                    first.append(m[0])
                    second.append(m[1])
            else:
                raise ValueError(f"Unexpected data type {type(y[0][0])}")

            scaling_factor = numpy.floor(numpy.log10((max(first)))) - numpy.floor(
                numpy.log10((max(second)))
            )
            second = [c * (10**scaling_factor) for c in second]
            points = {
                f"First Component - {self.mixture.first_component.name}": (
                    x,
                    first,
                    curve,
                ),
                f"Second Component - {self.mixture.second_component.name}, multiplied by {int(10**scaling_factor):.0e} ": (
                    x,
                    second,
                    curve,
                ),
            }

        elif type(y[0]) == float or type(y[0]) == numpy.float64:

            points = {f"{y_label}": (x, y, curve)}

        elif isinstance(y[0], Composition):

            points = {
                f"{y[0].type} fraction of {self.mixture.first_component.name}": (
                    x,
                    [c.first for c in y],
                    curve,
                )
            }
        else:
            raise ValueError(f"Unexpected data type {type(y[0])}")

        plot_graph(
            x_label=f"{self.mixture.first_component.name}, {self.feed_compositions[0].type} %",
            y_label=y_label,
            points=points,
        )

        return

    @property
    def get_separation_factor(self) -> typing.List[float]:
        """
        :return: List of separation Factors
        """
        feed = self.feed_compositions
        permeate = self.permeate_composition
        return [
            (permeate[i].first / permeate[i].second) / (feed[i].first / feed[i].second)
            for i in range(len(feed))
        ]

    @property
    def get_psi(self) -> typing.List[float]:
        """
        :return: List of Pervaporation Separation Index (PSI) values
        """
        separation_factor = self.get_separation_factor
        total_flux = [
            sum(self.partial_fluxes[i]) for i in range(len(self.partial_fluxes))
        ]
        return [numpy.multiply(total_flux, numpy.subtract(separation_factor, 1))]

    @property
    def get_selectivity(self) -> typing.List[float]:
        """
        :return: List of Calculated selectivities
        """
        permeance = self.permeances
        return [
            permeance[i][0].value / permeance[i][1].value for i in range(len(permeance))
        ]

    @classmethod
    def from_csv(cls, path: typing.Union[str, Path]) -> "ProcessModel":
        pass

    def to_csv(self, path: typing.Union[str, Path] = results_path) -> None:

        # mixture: Mixture
        # membrane_name: str +
        # feed_temperature: typing.List[float] +
        # feed_compositions: typing.List[Composition] ++
        # permeate_composition: typing.List[Composition] ++
        # permeate_temperature: typing.List[float] +
        # permeate_pressure: typing.List[float] +
        # feed_mass: typing.List[float] +
        # partial_fluxes: typing.List[typing.Tuple[float, float]] ++
        # permeances: typing.List[typing.Tuple[Permeance, Permeance]] +++
        # time: typing.List[float] +
        # feed_evaporation_heat: typing.List[float] +
        # permeate_condensation_heat: typing.List[typing.Optional[float]] +
        # initial_conditions: Conditions
        # permeance_fits: typing.Optional[
        #     typing.Tuple[PervaporationFunction, PervaporationFunction]
        # ] = None
        # comments: typing.Optional[str] = None

        process_dir = (
            path
            / f"ProcessModel_ID_{randrange(1,100)}_{datetime.now().minute}{datetime.now().second}"
        )

        process_dir.mkdir(parents=True, exist_ok=False)

        output = pandas.DataFrame(
            {
                "feed_temperature": [t for t in self.feed_temperature],
                "time": [t for t in self.time],
                "composition": [c.first for c in self.feed_compositions],
                "composition_type": [c.type for c in self.feed_compositions],
                "permeate_composition": [c.first for c in self.permeate_composition],
                "permeate_composition_type": [
                    c.type for c in self.permeate_composition
                ],
                "permeate_temperature": [p_t for p_t in self.permeate_temperature],
                "permeate_pressure": [p_p for p_p in self.permeate_pressure],
                "feed_mass": [m for m in self.feed_mass],
                "partial_flux_1": [f[0] for f in self.partial_fluxes],
                "partial_flux_2": [f[1] for f in self.partial_fluxes],
                "permeance_1": [p[0].value for p in self.permeances],
                "permeance_2": [p[1].value for p in self.permeances],
                "units": [p[0].units for p in self.permeances],
                "feed_evaporation_heat": [h for h in self.feed_evaporation_heat],
                "permeate condensation_heat": [
                    c_h for c_h in self.permeate_condensation_heat
                ],
            }
        )
        output["membrane_name"] = self.membrane_name
        output["mixture"] = self.mixture.name
        output["comment"] = self.comments

        output = output[PROCESS_MODEL_COLUMNS]

        output.to_csv(process_dir, index=False)

        self.permeance_fits[0].save((process_dir / f"PervaporationFunction_{self.mixture.first_component.name}.pv"))
        self.permeance_fits[1].save((process_dir / f"PervaporationFunction_{self.mixture.second_component.name}.pv"))
        joblib.dump(self.initial_conditions, (process_dir / "Initial_Conditions"))
        pass
