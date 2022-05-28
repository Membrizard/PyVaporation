import typing

import attr
import numpy

from diffusion_curve import DiffusionCurve
from process import ProcessModel
from conditions import Conditions
from membrane import Membrane
from mixture import Composition, CompositionType, Mixture, get_nrtl_partial_pressures


def non_ideal_model():
    return 0
