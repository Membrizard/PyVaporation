import attr

from enum import Enum

from ..mixture import Mixture

R = 8.314462


class CompositionType(Enum):
    molar: str = "molar"
    weight: str = "weight"


@attr.s(auto_attribs=True)
class Composition:
    p: float = attr.ib(validator=lambda value: 0 <= value <= 1)  # type: ignore
    type: CompositionType

    @property
    def first(self) -> float:
        return self.p

    @property
    def second(self) -> float:
        return 1 - self.p

    def to_molar(self, mixture: Mixture) -> "Composition":
        if self.type == CompositionType.molar:
            return self
        else:
            p = (self.p / mixture.first_component.molecular_weight) / (
                    self.p / mixture.first_component.molecular_weight
                    + (1 - self.p) / mixture.second_component.molecular_weight
            )
            return Composition(p=p, type=CompositionType('weight'))

    def to_weight(self, mixture: Mixture) -> "Composition":
        if self.type == CompositionType.weight:
            return self
        else:
            p = (mixture.first_component.molecular_weight * self.p) / (
                mixture.first_component.molecular_weight * self.p
                + mixture.second_component.molecular_weight * (1 - self.p)
            )
            return Composition(p=p, type=CompositionType('weight'))


@attr.s(auto_attribs=True)
class AntoineConstants:
    a: float = attr.ib(converter=lambda value: float(value))  # type: ignore
    b: float = attr.ib(converter=lambda value: float(value))  # type: ignore
    c: float = attr.ib(converter=lambda value: float(value))  # type: ignore


@attr.s(auto_attribs=True)
class NRTLParameters:
    g12: float
    g21: float
    alpha: float
