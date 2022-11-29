import typing
from datetime import datetime
from pathlib import Path

import attr
import joblib
import numpy
import pandas

from ..conditions import Conditions
from ..mixtures import Composition, Mixture, Mixtures
from ..optimizer import PervaporationFunction
from ..permeance import Permeance, Units
from ..plotting import plot_graph

PROCESS_MODEL_COLUMNS = [
    "membrane_name",
    "mixture",
    "time",
    "feed_mass",
    "feed_temperature",
    "permeate_temperature",
    "permeate_pressure",
    "composition",
    "composition_type",
    "permeate_composition",
    "permeate_composition_type",
    "partial_flux_1",
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
    """
    Class for description, working with and storage of process models
    """

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
    membrane_path: typing.Optional[Path] = None

    @staticmethod
    def _generate_process_path(membrane_path: typing.Union[str, Path]) -> Path:
        """
        :param membrane_path: path to an associated membrane
        :return: generates a path to save a process model
        """
        if type(membrane_path) is not Path:
            membrane_path = Path(membrane_path)

        results_path = membrane_path / "results"
        results_path.mkdir(parents=True, exist_ok=True)

        process_path = results_path / ("process_" + str(hash(datetime.now()))[0:4])
        process_path.mkdir(parents=True, exist_ok=False)

        return process_path

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
        return list(numpy.multiply(total_flux, numpy.subtract(separation_factor, 1)))

    @property
    def get_selectivity(self) -> typing.List[float]:
        """
        :return: List of Calculated selectivities
        """
        permeance = self.permeances
        return [
            permeance[i][0].value / permeance[i][1].value for i in range(len(permeance))
        ]

    # TODO Create safe_save  and save_load methods to save to and load from json objects

    @classmethod
    def load(
        cls,
        process_path: typing.Union[str, Path],
        is_safe: typing.Optional[bool] = False,
    ) -> "ProcessModel":
        """
        Creates a ProcessModel object from the data in a specified directory
        :param process_path: path to the process folder of a unified format
        :param is_safe: if True loads pv_func and initial condition objects from jsons
        :return: ProcessModel object
        """
        if type(process_path) is not Path:
            process_path = Path(process_path)

        process_frame = pandas.read_csv(process_path / "process_model.csv")

        pv_0_filenames = list(
            filter(
                lambda x: x.stem.startswith("pervaporation_function_0"),
                process_path.iterdir(),
            )
        )
        pv_1_filenames = list(
            filter(
                lambda x: x.stem.startswith("pervaporation_function_1"),
                process_path.iterdir(),
            )
        )
        if len(pv_0_filenames) != 1:
            raise FileExistsError(
                "Distinct first pervaporation function file is not found in %s"
                % process_path
            )
        if len(pv_1_filenames) != 1:
            raise FileExistsError(
                "Distinct second pervaporation function file is not found in %s"
                % process_path
            )

        if is_safe:
            pv_0 = PervaporationFunction.safe_load(pv_0_filenames[0])
            pv_1 = PervaporationFunction.safe_load(pv_1_filenames[0])
        else:
            pv_0 = PervaporationFunction.load(pv_0_filenames[0])
            pv_1 = PervaporationFunction.load(pv_1_filenames[0])

        if (process_path / "initial_conditions.ic").exists():
            if is_safe:
                initial_conditions = Conditions.safe_load(
                    process_path / "initial_conditions.ic"
                )
            else:
                initial_conditions = joblib.load(process_path / "initial_conditions.ic")
        else:
            initial_conditions = None

        mixture = getattr(Mixtures, process_frame["mixture"].iloc[0])

        if pandas.isna(process_frame["permeate_temperature"].iloc[0]):
            permeate_temperature = None
        else:
            permeate_temperature = process_frame["permeate_temperature"].iloc[0]
        if pandas.isna(process_frame["permeate_pressure"].iloc[0]):
            permeate_pressure = None
        else:
            permeate_pressure = process_frame["permeate_pressure"].iloc[0]

        if (
            process_frame["partial_flux_1"].isna().mean() == 0
            and process_frame["partial_flux_2"].isna().mean() == 0
        ):
            partial_fluxes = []
            for i in range(len(process_frame)):
                partial_fluxes.append(
                    (
                        process_frame["partial_flux_1"].iloc[i],
                        process_frame["partial_flux_2"].iloc[i],
                    )
                )
        else:
            partial_fluxes = None

        if (
            process_frame["permeance_1"].isna().mean() == 0
            and process_frame["permeance_2"].isna().mean() == 0
            and process_frame["units"].isna().mean() == 0
        ):
            permeances = []
            for i in range(len(process_frame)):
                permeances.append(
                    (
                        Permeance(
                            value=process_frame["permeance_1"].iloc[i],
                            units=process_frame["units"].iloc[0],
                        ).convert(
                            to_units=Units.kg_m2_h_kPa,
                            component=mixture.first_component,
                        ),
                        Permeance(
                            value=process_frame["permeance_2"].iloc[i],
                            units=process_frame["units"].iloc[0],
                        ).convert(
                            to_units=Units.kg_m2_h_kPa,
                            component=mixture.second_component,
                        ),
                    )
                )
        else:
            permeances = None

        return cls(
            mixture=mixture,
            membrane_name=process_frame["membrane_name"].iloc[0],
            feed_temperature=process_frame["feed_temperature"],
            feed_compositions=[
                Composition(
                    p=process_frame["composition"].iloc[i],
                    type=process_frame["composition_type"].iloc[i],
                ).to_weight(mixture=mixture)
                for i in range(len(process_frame))
            ],
            permeate_composition=[
                Composition(
                    p=process_frame["permeate_composition"].iloc[i],
                    type=process_frame["permeate_composition_type"].iloc[i],
                ).to_weight(mixture=mixture)
                for i in range(len(process_frame))
            ],
            permeate_temperature=permeate_temperature,
            permeate_pressure=permeate_pressure,
            feed_mass=process_frame["feed_mass"],
            partial_fluxes=partial_fluxes,
            permeances=permeances,
            time=process_frame["time"],
            feed_evaporation_heat=process_frame["feed_evaporation_heat"],
            permeate_condensation_heat=process_frame["permeate_condensation_heat"],
            initial_conditions=initial_conditions,
            permeance_fits=(pv_0, pv_1)
            if (pv_1 is not None and pv_0 is not None)
            else None,
            comments=process_frame["comment"],
            membrane_path=None,
        )

    def save(
        self,
        membrane_path: typing.Union[str, Path] = membrane_path,
        is_safe: typing.Optional[bool] = False,
    ) -> None:
        """
        Saves a ProcessModel object to a specified directory of a unified format
        :param membrane_path: path to the associated membrane
        :param is_safe: if True saves pv_func and initial condition objects to jsons
        :return: saves the object
        """
        if type(membrane_path) is not Path:
            membrane_path = Path(membrane_path)

        results_path = membrane_path / "results"
        results_path.mkdir(parents=True, exist_ok=True)

        process_path = self._generate_process_path(membrane_path)

        process_frame = pandas.DataFrame(
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
                "permeate_condensation_heat": [
                    c_h for c_h in self.permeate_condensation_heat
                ],
            }
        )
        process_frame["membrane_name"] = self.membrane_name
        process_frame["mixture"] = self.mixture.name
        process_frame["comment"] = self.comments
        process_frame = process_frame[PROCESS_MODEL_COLUMNS]
        process_frame.to_csv(process_path / "process_model.csv", index=False)

        if self.permeance_fits is None:
            self.permeance_fits = (PervaporationFunction(n=0, m=0, alpha=self.permeances[0][0].value, a=[0], b=[0]),
                                   PervaporationFunction(n=0, m=0, alpha=self.permeances[0][1].value, a=[0], b=[0]),)

        if is_safe:
            self.permeance_fits[0].safe_save(
                (
                    process_path
                    / f"pervaporation_function_0_{self.mixture.first_component.name}.pv"
                )
            )

            self.permeance_fits[1].safe_save(
                (
                    process_path
                    / f"pervaporation_function_1_{self.mixture.second_component.name}.pv"
                )
            )

            self.initial_conditions.safe_save(process_path / "initial_conditions.ic")
        else:

            self.permeance_fits[0].save(
                (
                    process_path
                    / f"pervaporation_function_0_{self.mixture.first_component.name}.pv"
                )
            )
            self.permeance_fits[1].save(
                (
                    process_path
                    / f"pervaporation_function_1_{self.mixture.second_component.name}.pv"
                )
            )

            joblib.dump(
                self.initial_conditions, (process_path / "initial_conditions.ic")
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
                f"Second Component - {self.mixture.second_component.name}, "
                f"multiplied by {int(10**scaling_factor):.0e} ": (
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
            x_label=f"Process time, hours",
            y_label=y_label,
            points=points,
        )
