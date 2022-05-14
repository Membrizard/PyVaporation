import typing

import attr
import numpy


# @attr.s(auto_attribs=True)
# class Membrane:
#     ideal_experiments: typing.List[IdealExperiment]
#     non_ideal_experiments: typing.List[NonIdealExperiment]
#
#     @classmethod
#     def from_ideal_csv(cls, path) -> 'Membrane':
#         pass
#
#     @classmethod
#     def from_non_ideal_csv(cls, path) -> 'Membrane':
#         pass
#
#     @property
#     def ideal_experiments_names(self) -> typing.List[str]:
#         return [ie.component.name for ie in self.ideal_experiments]
#
#     # Get all the penetrants the membrane was tested for
#     def get_known_penetrants(self) -> typing.List[Component]:
#         return numpy.unique(
#             [ideal_experiment.component for ideal_experiment in self.ideal_experiments]
#         )
#
#     # Picking only Experiments related to the chosen component
#     def get_penetrant_data(self, component) -> typing.List[IdealExperiment]:
#         return list(
#             filter(
#                 lambda value: value.component.name in self.ideal_experiments_names,
#                 self.ideal_experiments,
#             )
#         )
#
#     # Calculate an apparent activation energy of permeation
#     def calculate_activation_energy(self, component) -> float:
#         self.get_penetrant_data(component)
#         # Calculation of Least-squares Linear Fit of Ln(Permeance) versus 1/T
#         return 0
#
#     def get_permeance(self, temperature, component) -> float:
#         return 0
#
#     def get_ideal_selectivity(self, temperature, component1, component2) -> float:
#         return 0
