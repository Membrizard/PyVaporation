import attr

from mixture import Composition


@attr.s(auto_attribs=True)
class Conditions:
    membrane_area: float
    feed_temperature: float
    permeate_temperature: float
    feed_amount: float
    initial_feed_composition: Composition
