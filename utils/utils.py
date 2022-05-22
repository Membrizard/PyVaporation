import typing

import attr

R = 8.314462


@attr.s(auto_attribs=True)
class AntoineConstants:
    a: float = attr.ib(converter=lambda value: float(value))  # type: ignore
    b: float = attr.ib(converter=lambda value: float(value))  # type: ignore
    c: float = attr.ib(converter=lambda value: float(value))  # type: ignore


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
