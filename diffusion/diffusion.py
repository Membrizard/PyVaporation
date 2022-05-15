@attr.s(auto_attribs=True)
class DiffusionCurve:
    feed_temperature: float
    permeate_temperature: typing.Optional[float] = 0
    mixture: Mixture
    composition_range: typing.List[Composition]
    partial_fluxes: typing.List[typing.List[float]]
    total_flux: typing.List[float]
    separation_factor: typing.List[float]
    PSI: typing.List[float]
