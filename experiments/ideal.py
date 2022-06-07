import typing
from pathlib import Path

import attr
import pandas

from component import AllComponents, Component
from permeance import Permeance

# def convert_component(name: str, all_components: AllComponents) -> Component:
#     return getattr(all_components, name)


@attr.s(auto_attribs=True)
class IdealExperiment:
    """
    Class for specification of Ideal experiments, where Permeance is assumed to be constant
    over considered composition range
    """

    name: str
    temperature: float
    component: Component  # TODO: use real components
    permeance: Permeance
    activation_energy: typing.Optional[float] = None
    comment: typing.Optional[str] = None

    @classmethod
    def from_dict(
        cls,
        d: typing.Mapping[str, typing.Union[str, float]],
        all_components: AllComponents,
    ) -> "IdealExperiment":
        component = getattr(all_components, d["component"])
        return cls(**d)

    # TODO Add check for units and conversion to kg/(m2*h*kPa) using Permeance.convert()


@attr.s(auto_attribs=True)
class IdealExperiments:
    experiments: typing.List[IdealExperiment]

    def __len__(self):
        return len(self.experiments)

    @classmethod
    def from_csv(cls, path: typing.Union[str, Path]) -> "IdealExperiments":
        frame = pandas.read_csv(path)

        experiments = []
        for _, row in frame.iterrows():
            experiments.append(IdealExperiment.from_dict(row.to_dict()))

        return IdealExperiments(experiments=experiments)
