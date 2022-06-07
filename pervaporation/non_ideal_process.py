import typing

import attr
import numpy

from conditions import Conditions
from diffusion_curve import DiffusionCurve
from membrane import Membrane
from mixture import Composition, CompositionType, Mixture, get_nrtl_partial_pressures
from process import ProcessModel


def non_ideal_model():
    return 0
