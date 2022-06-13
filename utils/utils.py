import typing

import attr

R = 8.314462


class VPConstantsType:
    antoine: str = "antoine"
    frost: str = "frost"


@attr.s(auto_attribs=True)
class VaporPressureConstants:
    a: float = attr.ib(converter=lambda value: float(value))  # type: ignore
    b: float = attr.ib(converter=lambda value: float(value))  # type: ignore
    c: float = attr.ib(converter=lambda value: float(value))  # type: ignore
    type: str = attr.ib(
        default=VPConstantsType.antoine,
        converter=lambda x: getattr(VPConstantsType, x)
        if x is not None
        else VPConstantsType.antoine,
    )  # type: ignore


@attr.s(auto_attribs=True)
class NRTLParameters:
    g12: float
    g21: float
    alpha12: float
    alpha21: typing.Optional[float] = None
    a12: typing.Optional[float] = 0
    a21: typing.Optional[float] = 0


@attr.s(auto_attribs=True)
class HeatCapacityConstants:
    a: float = attr.ib(converter=lambda value: float(value))  # type: ignore
    b: float = attr.ib(converter=lambda value: float(value))  # type: ignore
    c: float = attr.ib(converter=lambda value: float(value))  # type: ignore
    d: float = attr.ib(converter=lambda value: float(value))  # type: ignore
