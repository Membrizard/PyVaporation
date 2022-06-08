import typing

import attr
import numpy

from conditions import Conditions
from diffusion_curve import DiffusionCurve
from membrane import Membrane
from mixture import Composition, CompositionType, Mixture, get_nrtl_partial_pressures
from process import ProcessModel
from optimizer import PervaporationFunction


def non_ideal_model(
        permeance_function_first: PervaporationFunction,
        permeance_function_second: PervaporationFunction,

):
    return 0
