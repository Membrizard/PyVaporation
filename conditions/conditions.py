import typing

from mixture import Composition

import attr


@attr.s(auto_attribs=True)
class Conditions:
    membrane_area: float
    feed_temperature: float
    feed_amount: float
    initial_feed_composition: Composition
    permeate_temperature: float
