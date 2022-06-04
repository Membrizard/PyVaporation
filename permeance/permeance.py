import attr
import typing
from component import Component


@attr.s(auto_attribs=True)
class Units:
    GPU: str = "GPU"
    SI: str = "SI"
    kg_m2_h_kPa: str = "kg/(m2*h*kPa)"
    Barrer: str = "Barrer"


@attr.s(auto_attribs=True)
class Permeance:
    value: float
    units: str = Units().kg_m2_h_kPa

    def convert(
        self,
        to_units: str,
        component: typing.Optional[Component] = None,
    ) -> "Permeance":

        if component is None and to_units == Units().kg_m2_h_kPa:
            raise ValueError("For conversion from and to kg/(m2*h*kPa) Component must be specified ")
        try:
            if component is None:
                conversion_dict = {
                    "GPU": 3.35e-10,
                    "SI": 1,
                    "Barrer": 3.35e-16,
                }
            else:
                conversion_dict = {
                    "GPU": 3.35e-10,
                    "SI": 1,
                    "Barrer": 3.35e-16,
                    "kg/(m2*h*kPa)": component.molecular_weight * 3.6e3,
                }
        except KeyError:
            raise KeyError("Conversion to specified units is not supported")

        return Permeance(
                value=(
                        self.value * conversion_dict[self.units] / conversion_dict[to_units]
                ),
                units=to_units,
            )

