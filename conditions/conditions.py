import attr

from mixture import Composition


@attr.s(auto_attribs=True)
class Conditions:
    membrane_area: float
    temperature: float
    permeate_temperature: float
    feed_amount: float
    feed_composition: Composition
    isothermal: bool = True
