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
class UNIQUACConstants:
    r: float
    q_geometric: float
    q_interaction: typing.Optional[float]

    def __attrs_post_init__(self):
        if self.q_interaction is None:
            self.q_interaction = self.q_geometric


@attr.s(auto_attribs=True)
class UNIQUACParameters:
    alpha_12: float
    alpha_21: float
    beta_12: float
    beta_21: float
    z: int = 10


@attr.s(auto_attribs=True)
class HeatCapacityConstants:
    a: float = attr.ib(converter=lambda value: float(value))  # type: ignore
    b: float = attr.ib(converter=lambda value: float(value))  # type: ignore
    c: float = attr.ib(converter=lambda value: float(value))  # type: ignore
    d: float = attr.ib(converter=lambda value: float(value))  # type: ignore
