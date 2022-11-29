import typing
from pathlib import Path

import attr
import numpy
import pandas

from ..mixtures import (
    Composition,
    CompositionType,
    Mixture,
    Mixtures,
    get_nrtl_partial_pressures,
)
from ..permeance import Permeance, Units
from ..plotting import plot_graph

DC_SET_COLUMNS = [
    "curve_id",
    "membrane_name",
    "mixture",
    "feed_temperature",
    "permeate_temperature",
    "permeate_pressure",
    "composition",
    "composition_type",
    "partial_flux_1",
    "partial_flux_2",
    "permeance_1",
    "permeance_2",
    "units",
    "comment",
]


@attr.s(auto_attribs=True)
class DiffusionCurve:
    """
    The class for Diffusion Curves - dependency of selective-transport properties
    on feed composition, and process temperature
    """

    mixture: Mixture
    membrane_name: str
    feed_temperature: float
    feed_compositions: typing.List[Composition]
    partial_fluxes: typing.Optional[typing.List[typing.Tuple[float, float]]] = None
    permeate_temperature: typing.Optional[float] = None
    permeate_pressure: typing.Optional[float] = None
    permeances: typing.Optional[typing.List[typing.Tuple[Permeance, Permeance]]] = None
    comments: typing.Optional[str] = None

    def __attrs_post_init__(self):

        if self.partial_fluxes is None and self.permeances is not None:
            """
            Calculation of Partial Fluxes values if they are not specified
            Permeance values are converted to kg/(m2*h*kPa),
            The permeate temperature and pressure parameters are ignored, test_components' pressure in permeate
            is considered zero
            :return a list of Partial fluxes for each component tuple(Ji,Jj) at each concentration
            """
            feed_partial_pressures = [
                get_nrtl_partial_pressures(
                    self.feed_temperature, self.mixture, composition
                )
                for composition in self.feed_compositions
            ]

            self.permeances = [
                (
                    self.permeances[i][0].convert(
                        to_units=Units.kg_m2_h_kPa,
                        component=self.mixture.first_component,
                    ),
                    self.permeances[i][1].convert(
                        to_units=Units.kg_m2_h_kPa,
                        component=self.mixture.second_component,
                    ),
                )
                for i in range(len(self.permeances))
            ]

            self.partial_fluxes = [
                (
                    self.permeances[i][0].value * feed_partial_pressures[i][0],
                    self.permeances[i][1].value * feed_partial_pressures[i][1],
                )
                for i in range(len(self.permeances))
            ]

        if self.permeances is None and self.partial_fluxes is not None:
            """
            Calculation of Permeance values if they are not specified
            Permeance values are in kg/(m2*h*kPa)
            If needed, Permeance values may be converted using Permeance.converter()
            :return a list of Permeances for each component tuple(Pi,Pj) at each concentration
            """
            permeate_compositions = self.permeate_composition
            feed_partial_pressures = [
                get_nrtl_partial_pressures(
                    self.feed_temperature, self.mixture, composition
                )
                for composition in self.feed_compositions
            ]
            if self.permeate_temperature is None and self.permeate_pressure is None:
                self.permeances = [
                    (
                        Permeance(
                            value=self.partial_fluxes[i][0]
                            / feed_partial_pressures[i][0]
                        ),
                        Permeance(
                            value=self.partial_fluxes[i][1]
                            / feed_partial_pressures[i][1]
                        ),
                    )
                    for i in range(len(self.feed_compositions))
                ]
            elif (
                self.permeate_temperature is not None and self.permeate_pressure is None
            ):
                permeate_partial_pressures = [
                    get_nrtl_partial_pressures(
                        self.permeate_temperature, self.mixture, composition
                    )
                    for composition in permeate_compositions
                ]
                self.permeances = [
                    (
                        Permeance(
                            value=self.partial_fluxes[i][0]
                            / (
                                feed_partial_pressures[i][0]
                                - permeate_partial_pressures[i][0]
                            )
                        ),
                        Permeance(
                            value=self.partial_fluxes[i][1]
                            / (
                                feed_partial_pressures[i][1]
                                - permeate_partial_pressures[i][1]
                            )
                        ),
                    )
                    for i in range(len(self.feed_compositions))
                ]

            elif (
                self.permeate_pressure is not None and self.permeate_temperature is None
            ):
                permeate_partial_pressures = [
                    (
                        self.permeate_pressure
                        * composition.to_molar(self.mixture).first,
                        self.permeate_pressure
                        * composition.to_molar(self.mixture).second,
                    )
                    for composition in permeate_compositions
                ]
                self.permeances = [
                    (
                        Permeance(
                            value=self.partial_fluxes[i][0]
                            / (
                                feed_partial_pressures[i][0]
                                - permeate_partial_pressures[i][0]
                            )
                        ),
                        Permeance(
                            value=self.partial_fluxes[i][1]
                            / (
                                feed_partial_pressures[i][1]
                                - permeate_partial_pressures[i][1]
                            )
                        ),
                    )
                    for i in range(len(self.feed_compositions))
                ]
            else:
                raise ValueError(
                    "Either permeate temperature or permeate pressure could be stated not both"
                )
        elif self.permeances is None and self.partial_fluxes is None:
            raise ValueError(
                "Either Permeances or Fluxes must be specified as functions of feed composition"
            )
        else:
            self.permeances = [
                (
                    self.permeances[i][0].convert(
                        to_units=Units.kg_m2_h_kPa,
                        component=self.mixture.first_component,
                    ),
                    self.permeances[i][1].convert(
                        to_units=Units.kg_m2_h_kPa,
                        component=self.mixture.second_component,
                    ),
                )
                for i in range(len(self.permeances))
            ]

    def __len__(self):
        return len(self.feed_compositions)

    def plot(self, y: typing.List, y_label: str = "", curve: bool = 1):
        """
        Draws basic plot of a specified parameter versus first component fraction in the feed
        :param y: parameter
        :param y_label: Name of the axis
        :param curve: if True - connects points with a line, if False - draws raw points
        :return: plots a graph
        """

        x = [c.first for c in self.feed_compositions]
        if type(y[0]) == tuple:
            first = []
            second = []
            if isinstance(y[0][0], Permeance):
                for m in y:
                    first.append(m[0].value)
                    second.append(m[1].value)
            elif type(y[0][0]) == float or type(y[0][0]) == numpy.float64:
                for m in y:
                    first.append(m[0])
                    second.append(m[1])
            else:
                raise ValueError(f"Unexpected default_membranes type {type(y[0][0])}")

            scaling_factor = numpy.floor(numpy.log10((max(first)))) - numpy.floor(
                numpy.log10((max(second)))
            )
            second = [c * (10**scaling_factor) for c in second]
            points = {
                f"First Component - {self.mixture.first_component.name}": (
                    x,
                    first,
                    curve,
                ),
                f"Second Component - {self.mixture.second_component.name}, "
                f"multiplied by {int(10**scaling_factor):.0e} ": (
                    x,
                    second,
                    curve,
                ),
            }

        elif type(y[0]) == float or type(y[0]) == numpy.float64:

            points = {f"{y_label}": (x, y, curve)}

        elif isinstance(y[0], Composition):
            points = {
                f"{y[0].type} fraction of {self.mixture.first_component.name}": (
                    x,
                    [c.first for c in y],
                    curve,
                )
            }
        else:
            raise ValueError(f"Unexpected default_membranes type {type(y[0])}")

        plot_graph(
            x_label=f"{self.mixture.first_component.name}, {self.feed_compositions[0].type} %",
            y_label=y_label,
            points=points,
        )

    @property
    def permeate_composition(self) -> typing.List[Composition]:
        """
        Calculation of permeate compositions at each concentration
        :return - List of permeate weight compositions
        """
        return [
            Composition(
                self.partial_fluxes[i][0] / (sum(self.partial_fluxes[i])),
                type=CompositionType.weight,
            )
            for i in range(len(self.feed_compositions))
        ]

    @property
    def get_separation_factor(self) -> typing.List[float]:
        """
        Calculation of separation factor at each concentration
        :return - List of separation factors
        """
        permeate_composition = self.permeate_composition
        feed_composition = self.feed_compositions
        return [
            (permeate_composition[i].first / permeate_composition[i].second)
            / (feed_composition[i].first / feed_composition[i].second)
            for i in range(len(self.feed_compositions))
        ]

    @property
    def get_psi(self) -> typing.List[float]:
        """
        Calculation of Pervaporation Separation Index (PSI) at each concentration
        :return - List of PSI
        """
        total_flux = [
            sum(self.partial_fluxes[i]) for i in range(len(self.feed_compositions))
        ]
        separation_factor = self.get_separation_factor
        return list(numpy.multiply(total_flux, numpy.subtract(separation_factor, 1)))

    @property
    def get_permeances(self) -> typing.List[typing.Tuple[Permeance, Permeance]]:
        """
        Calculates permenace of both components based on the partial flux values, mixture, feed  and permeate parameters
        :return: Permeances of each component at each feed concentration as [(P1,P2),...]
        """
        if self.permeances is not None:
            return self.permeances
        else:
            permeate_compositions = self.permeate_composition
            feed_partial_pressures = [
                get_nrtl_partial_pressures(
                    self.feed_temperature, self.mixture, composition
                )
                for composition in self.feed_compositions
            ]
            if self.permeate_temperature is None:
                self.permeances = [
                    (
                        Permeance(
                            value=self.partial_fluxes[i][0]
                            / feed_partial_pressures[i][0]
                        ),
                        Permeance(
                            value=self.partial_fluxes[i][1]
                            / feed_partial_pressures[i][1]
                        ),
                    )
                    for i in range(len(self.feed_compositions))
                ]
            else:
                permeate_partial_pressures = [
                    get_nrtl_partial_pressures(
                        self.permeate_temperature, self.mixture, composition
                    )
                    for composition in permeate_compositions
                ]
                self.permeances = [
                    (
                        Permeance(
                            value=self.partial_fluxes[i][0]
                            / (
                                feed_partial_pressures[i][0]
                                - permeate_partial_pressures[i][0]
                            )
                        ),
                        Permeance(
                            value=self.partial_fluxes[i][1]
                            / (
                                feed_partial_pressures[i][1]
                                - permeate_partial_pressures[i][1]
                            )
                        ),
                    )
                    for i in range(len(self.feed_compositions))
                ]

    @property
    def get_selectivity(
        self,
    ) -> typing.List[float]:
        """
        Calculation of selectivity at each concentration
        :return - List of selectivities (Permeances are in SI by default)
        """
        permeances = self.permeances
        return [
            permeances[i][0].convert("SI", self.mixture.first_component).value
            / permeances[i][1].convert("SI", self.mixture.second_component).value
            for i in range(len(self.feed_compositions))
        ]

    @classmethod
    def from_frame(cls, data: pandas.DataFrame) -> "DiffusionCurve":
        # TODO: docstring
        mixture = getattr(Mixtures, data["mixture"].iloc[0])

        if (
            data["partial_flux_1"].isna().mean() == 0
            and data["partial_flux_2"].isna().mean() == 0
        ):
            partial_fluxes = []
            for i in range(len(data)):
                partial_fluxes.append(
                    (
                        data["partial_flux_1"].iloc[i],
                        data["partial_flux_2"].iloc[i],
                    )
                )
        else:
            partial_fluxes = None

        if (
            data["permeance_1"].isna().mean() == 0
            and data["permeance_2"].isna().mean() == 0
            and data["units"].isna().mean() == 0
        ):
            permeances = []
            for i in range(len(data)):
                permeances.append(
                    (
                        Permeance(
                            value=data["permeance_1"].iloc[i],
                            units=data["units"].iloc[0],
                        ).convert(
                            to_units=Units.kg_m2_h_kPa,
                            component=mixture.first_component,
                        ),
                        Permeance(
                            value=data["permeance_2"].iloc[i],
                            units=data["units"].iloc[0],
                        ).convert(
                            to_units=Units.kg_m2_h_kPa,
                            component=mixture.second_component,
                        ),
                    )
                )
        else:
            permeances = None

        if partial_fluxes is None and permeances is None:
            raise ValueError("Partial fluxes and permeances are not provided")

        if pandas.isna(data["permeate_temperature"].iloc[0]):
            permeate_temperature = None
        else:
            permeate_temperature = data["permeate_temperature"].iloc[0]
        if pandas.isna(data["permeate_pressure"].iloc[0]):
            permeate_pressure = None
        else:
            permeate_pressure = data["permeate_pressure"].iloc[0]

        return DiffusionCurve(
            mixture=mixture,
            membrane_name=data["membrane_name"].iloc[0],
            feed_temperature=data["feed_temperature"].iloc[0],
            feed_compositions=[
                Composition(
                    p=data["composition"].iloc[i], type=data["composition_type"].iloc[i]
                ).to_weight(mixture=mixture)
                for i in range(len(data))
            ],
            partial_fluxes=partial_fluxes,
            permeate_temperature=permeate_temperature,
            permeate_pressure=permeate_pressure,
            permeances=permeances,
            comments=data["comment"].iloc[0],
        )

    def save(self, path: typing.Union[str, Path]) -> None:
        """
        Saves the curve to a specified path in the form of a .csv file
        :param path: path to the directory
        :return: None, the object is saved in a specified directory
        """
        output = pandas.DataFrame(
            {
                "composition": [c.p for c in self.feed_compositions],
                "composition_type": [c.type for c in self.feed_compositions],
                "partial_flux_1": [f[0] for f in self.partial_fluxes],
                "partial_flux_2": [f[1] for f in self.partial_fluxes],
                "permeance_1": [p[0].value for p in self.permeances],
                "permeance_2": [p[1].value for p in self.permeances],
                "units": [p[0].units for p in self.permeances],
            }
        )
        output["curve_id"] = "1"
        output["membrane_name"] = self.membrane_name
        output["mixture"] = self.mixture.name
        output["feed_temperature"] = self.feed_temperature
        output["permeate_temperature"] = self.permeate_temperature
        output["permeate_pressure"] = self.permeate_pressure
        output["comment"] = self.comments
        output = output[DC_SET_COLUMNS]
        output.to_csv(path, index=False)


@attr.s(auto_attribs=True)
class DiffusionCurveSet:
    """
    A class for storing and working with multiple diffusion curves related to the same set of experiments.
    """

    name: str
    diffusion_curves: typing.List[DiffusionCurve]

    def __getitem__(self, item):
        return self.diffusion_curves[item]

    @classmethod
    def load(cls, path: Path) -> "DiffusionCurveSet":
        data = pandas.read_csv(path)
        name = path.stem

        if list(data.columns) != DC_SET_COLUMNS:
            raise ValueError(
                "Incorrect default_membranes: %s at %s" % (list(data.columns), path)
            )

        raw_diffusion_curves = data.groupby("curve_id")
        diffusion_curves = []
        for _, curve in raw_diffusion_curves:
            diffusion_curves.append(DiffusionCurve.from_frame(curve))

        return cls(
            name=name,
            diffusion_curves=diffusion_curves,
        )
