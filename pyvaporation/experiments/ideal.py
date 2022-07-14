import typing
from pathlib import Path

import attr
import pandas

from ..components import Component, Components
from ..permeance import Permeance, Units

IDEAL_EXPERIMENT_COLUMNS = [
    "name",
    "temperature",
    "component",
    "activation_energy",
    "permeance",
    "units",
    "comment",
]


@attr.s(auto_attribs=True)
class IdealExperiment:
    """
    Class for specification of Ideal experiments, where Permeance is assumed to be constant
    over considered composition range
    """

    name: str
    temperature: float
    component: Component
    permeance: Permeance
    activation_energy: typing.Optional[float] = None
    comment: typing.Optional[str] = None

    @classmethod
    def from_dict(
        cls, d: typing.Mapping[str, typing.Union[str, float]]
    ) -> "IdealExperiment":
        """
        Reads IdealExperiment objects from dictionary
        :param d: Mapping of IdealExperiment parameters
        :return: IdealExperiment object with parameters specified in mapping
        """
        component = getattr(Components, d["component"])
        permeance = Permeance(value=d["permeance"], units=d["units"]).convert(
            to_units=Units.kg_m2_h_kPa, component=component
        )
        return cls(
            name=d["name"],
            temperature=d["temperature"],
            component=component,
            permeance=permeance,
            activation_energy=d["activation_energy"],
            comment=d["comment"],
        )


@attr.s(auto_attribs=True)
class IdealExperiments:
    """
    A class for storing multiple IdealExperiments related to a particular membrane
    """

    experiments: typing.List[IdealExperiment]

    def __len__(self):
        return len(self.experiments)

    @classmethod
    def from_csv(cls, path: typing.Union[str, Path]) -> "IdealExperiments":
        """
        Generates IdealExperiments object from the .csv file
        :param path: Path to the file
        :return: IdealExperiments object with IdealExperiments described in a .csv file
        """
        frame = pandas.read_csv(path)
        if list(frame.columns) != IDEAL_EXPERIMENT_COLUMNS:
            raise ValueError("Incorrect columns: %s" % list(frame.columns))

        experiments = []
        for _, row in frame.iterrows():
            experiments.append(IdealExperiment.from_dict(row.to_dict()))

        return IdealExperiments(experiments=experiments)
