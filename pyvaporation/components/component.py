import attr
import numpy

from ..utils import HeatCapacityConstants, R, VaporPressureConstants, VPConstantsType


@attr.s(auto_attribs=True)
class Component:
    name: str
    molecular_weight: float = attr.ib(converter=lambda value: float(value))
    vapour_pressure_constants: VaporPressureConstants
    heat_capacity_constants: HeatCapacityConstants

    def get_vapor_pressure(self, temperature: float) -> float:
        """
        Calculation of saturated pressure in kPa at a given temperature in K using Antoine equation (by the basis of 10)
        type of the calculation using Antoine (antoine) log10(P)=a+b/(T+C)
        or Frost equation (Frost) ln(P) = a+b/T+c/T^2
        :param temperature: temperature in K
        :return: saturated pressure in kPa calculated with respect to constants and given temperature
        """
        if self.vapour_pressure_constants.type == VPConstantsType.antoine:
            return 10 ** (
                self.vapour_pressure_constants.a
                + self.vapour_pressure_constants.b
                / (temperature + self.vapour_pressure_constants.c)
            )
        elif self.vapour_pressure_constants.type == VPConstantsType.frost:
            return numpy.exp(
                self.vapour_pressure_constants.a
                + self.vapour_pressure_constants.b / temperature
                + self.vapour_pressure_constants.c / temperature**2
            )
        else:
            raise ValueError("Type of calculation not supported")

    def get_vaporisation_heat(self, temperature: float) -> float:
        """
        Calculation of Vaporisation heat in kJ/mol using Clapeyron-Clausius equation
        for Antoine equation: H=b*R*ln(10)*(T/(T+C))^2
        for Frost Equation: H=R*(-b-2*c/T)
        :param temperature: temperature in K
        :return: Vaporisation heat in kJ/mol
        """
        if self.vapour_pressure_constants.type == VPConstantsType.antoine:
            return (
                -(
                    (temperature / (temperature + self.vapour_pressure_constants.c))
                    ** 2
                    * R
                    * self.vapour_pressure_constants.b
                    * numpy.log(10)
                )
                / 1000
            )
        elif self.vapour_pressure_constants.type == VPConstantsType.frost:
            return (
                -R
                * (
                    self.vapour_pressure_constants.b
                    + 2 * self.vapour_pressure_constants.c / temperature
                )
                / 1000
            )
        else:
            raise ValueError("Type of calculation not supported")

    def get_specific_heat(self, temperature: float) -> float:
        """
        Calculation of Heat Capacity in J/(mol*K) using polynomial isobaric heat capacity fit
        :param temperature: temperature in K
        :return: Isobaric Heat capacity in J/(mol*K)
        """
        return (
            self.heat_capacity_constants.a
            + self.heat_capacity_constants.b * temperature
            + self.heat_capacity_constants.c * temperature**2
            + self.heat_capacity_constants.d * temperature**3
        )

    def get_cooling_heat(self, t0, t1):
        """
        Calculation of Specific Heat in J/mol using Integral (T2-T1) (CpdT)
        :param t0: temperature in K (t1 < t0)
        :param t1: temperature in K (t1 < t0)
        :return: Specific Heat in J/mol
        """
        return (
            self.heat_capacity_constants.a * (t0 - t1)
            + self.heat_capacity_constants.b * (t0**2 - t1**2) / 2
            + self.heat_capacity_constants.c * (t0**3 - t1**3) / 3
            + self.heat_capacity_constants.d * (t0**4 - t1**4) / 4
        )
