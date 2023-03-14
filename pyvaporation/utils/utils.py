import typing
import numpy
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
    q_interaction: typing.Optional[float] = None

    def __attrs_post_init__(self):
        if self.q_interaction is None:
            self.q_interaction = self.q_geometric


@attr.s(auto_attribs=True)
class UNIQUACBinaryInteractionParameters:
    i_component: str
    j_component: str
    ij_parameter: typing.Tuple[float, float]
    ji_parameter: typing.Tuple[float, float]


@attr.s(auto_attribs=True)
class UNIQUACParameters:
    binary_parameters: typing.List[UNIQUACBinaryInteractionParameters]
    # alpha_12: float
    # alpha_21: float
    # beta_12: float
    # beta_21: float
    z: int = 10

    def __len__(self):
        return len(self.binary_parameters)

    @classmethod
    def from_array(
        cls, array: typing.Union[typing.List[float], numpy.ndarray]
    ) -> "UNIQUACParameters":
        assert len(array) == 5
        return cls(
            binary_parameters_matrix=[
                [0, (array[0], array[2])],
                [(array[1], array[3]), 0],
            ],
            # alpha_12=array[0],
            # alpha_21=array[1],
            # beta_12=array[2],
            # beta_21=array[3],
            z=int(array[4]),
        )

# @attr.s(auto_attribs=True)
# class UNIQUACParameters:
#     binary_parameters_matrix: typing.List[typing.List[typing.Tuple[float, float]]]
#     # alpha_12: float
#     # alpha_21: float
#     # beta_12: float
#     # beta_21: float
#     z: int = 10
#
#     def __len__(self):
#         return len(self.binary_parameters_matrix)
#
#     @classmethod
#     def from_array(
#         cls, array: typing.Union[typing.List[float], numpy.ndarray]
#     ) -> "UNIQUACParameters":
#         assert len(array) == 5
#         return cls(
#             binary_parameters_matrix=[
#                 [0, (array[0], array[2])],
#                 [(array[1], array[3]), 0],
#             ],
#             # alpha_12=array[0],
#             # alpha_21=array[1],
#             # beta_12=array[2],
#             # beta_21=array[3],
#             z=int(array[4]),
#         )


@attr.s(auto_attribs=True)
class HeatCapacityConstants:
    a: float = attr.ib(converter=lambda value: float(value))  # type: ignore
    b: float = attr.ib(converter=lambda value: float(value))  # type: ignore
    c: float = attr.ib(converter=lambda value: float(value))  # type: ignore
    d: float = attr.ib(converter=lambda value: float(value))  # type: ignore
