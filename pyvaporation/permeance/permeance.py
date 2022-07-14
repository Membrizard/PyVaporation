import typing

import attr

from ..components import Component


class Units:
    """
    Class to represent Permeance units
    """

    GPU: str = "GPU"
    SI: str = "SI"
    kg_m2_h_kPa: str = "kg/(m2*h*kPa)"


@attr.s(auto_attribs=True)
class Permeance:
    """
    Class to represent Permeance values
    """

    value: float = attr.ib(converter=lambda x: x if x >= 0 else 0)
    units: str = Units.kg_m2_h_kPa

    def __add__(self, other: "Permeance") -> "Permeance":
        if self.units != other.units:
            raise ValueError("Only Permeances in same units could be added")
        return Permeance(value=self.value + other.value, units=self.units)

    def convert(
        self,
        to_units: str,
        component: typing.Optional[Component] = None,
    ) -> "Permeance":
        """
        Converts Permeance values from and to commonly used units:
        kg/(m2*h*kPa),
        SI - (mol/(m2*s*Pa))
        GPU
        """
        if to_units == self.units:
            return self

        if component is None and to_units == Units().kg_m2_h_kPa:
            raise ValueError(
                "For conversion from and to kg/(m2*h*kPa) Component must be specified "
            )
        try:
            if component is None:
                conversion_dict = {
                    "GPU": 3.35e-10,
                    "SI": 1,
                }
            else:
                conversion_dict = {
                    "GPU": 3.35e-10,
                    "SI": 1,
                    "kg/(m2*h*kPa)": 1 / (component.molecular_weight * 3.6e3),
                }
        except KeyError:
            raise KeyError("Conversion to specified units is not supported")

        return Permeance(
            value=(
                self.value * conversion_dict[self.units] / conversion_dict[to_units]
            ),
            units=to_units,
        )
