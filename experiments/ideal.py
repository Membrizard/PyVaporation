import attr
import typing
import pandas

from pathlib import Path

from component import Component, AllComponents


# def convert_component(name: str, all_components: AllComponents) -> Component:
#     return getattr(all_components, name)


@attr.s(auto_attribs=True)
class IdealExperiment:
    # temperature: float  # K
    # permeance: float  # kg * mcm / (m2 * h * kPa)
    # component: Component
    # activation_energy: typing.Optional[float] = None

    name: str
    temperature: float
    component: str
    activation_energy: float
    permeance: float
    comment: str

    @classmethod
    def from_dict(cls, d: typing.Mapping[str, typing.Union[str, float]]) -> 'IdealExperiment':
        return cls(**d)


@attr.s(auto_attribs=True)
class IdealExperiments:
    experiments: typing.List[IdealExperiment]

    @classmethod
    def from_csv(cls, path: typing.Union[str, Path]) -> 'IdealExperiments':
        frame = pandas.read_csv(path)

        experiments = []
        for _, row in frame.iterrows():
            experiments.append(IdealExperiment.from_dict(row.to_dict()))

        return IdealExperiments(experiments=experiments)
