from pyvaporation.utils import HeatCapacityConstants, VaporPressureConstants

from .component import Component


class Components:
    """
    A class with pre-defined components with constants validated against literature data
    For details see ./tests/test_components/test_default_components.py
    """

    H2O: Component = Component(
        name="H2O",
        molecular_weight=18.02,
        vapour_pressure_constants=VaporPressureConstants(
            a=7.20389,
            b=-1733.926,
            c=-39.485,
        ),
        heat_capacity_constants=HeatCapacityConstants(
            a=32.2,
            b=1.924e-3,
            c=1.055e-5,
            d=-3.596e-9,
        ),
    )

    MeOH: Component = Component(
        name="MeOH",
        molecular_weight=32.04,
        vapour_pressure_constants=VaporPressureConstants(
            a=7.2209903,
            b=-1590.15535,
            c=-32.77001,
        ),
        heat_capacity_constants=HeatCapacityConstants(
            a=89.26344,
            b=-1.66291e-1,
            c=1.84499e-4,
            d=9.44251e-7,
        ),
    )

    EtOH: Component = Component(
        name="EtOH",
        molecular_weight=46.07,
        vapour_pressure_constants=VaporPressureConstants(
            a=7.24677,
            b=-1598.673,
            c=-46.424,
        ),
        heat_capacity_constants=HeatCapacityConstants(
            a=147.815652,
            b=-0.6732612305,
            c=0.001889017424,
            d=0,
        ),
    )

    iPOH: Component = Component(
        name="iPOH",
        molecular_weight=60.10,
        vapour_pressure_constants=VaporPressureConstants(
            a=6.861,
            b=-1357.427,
            c=-75.814,
        ),
        heat_capacity_constants=HeatCapacityConstants(
            a=53.8206447,
            b=7.92412193e-1,
            c=-4.56017218e-3,
            d=1.01887619e-5,
        ),
    )

    MTBE: Component = Component(
        name="MTBE",
        molecular_weight=88.15,
        vapour_pressure_constants=VaporPressureConstants(
            a=6.050931522,
            b=-1139.816725,
            c=-46.15171,
        ),
        heat_capacity_constants=HeatCapacityConstants(
            a=147.329712,
            b=-9.7850807e-2,
            c=9.216480e-4,
            d=-4.75200e-7,
        ),
    )

    ETBE: Component = Component(
        name="ETBE",
        molecular_weight=102.17,
        vapour_pressure_constants=VaporPressureConstants(
            a=6.0703,
            b=-1206.874,
            c=-49.19,
        ),
        heat_capacity_constants=HeatCapacityConstants(
            a=18.79,
            b=1.251,
            c=-3.015e-3,
            d=3.637e-6,
        ),
    )

    DME: Component = Component(
        name="DME",
        molecular_weight=90.12,
        vapour_pressure_constants=VaporPressureConstants(
            a=5.83775,
            b=-1260.52,
            c=-37.322,
        ),
        heat_capacity_constants=HeatCapacityConstants(
            a=1.43841e2,
            b=3.88448e-1,
            c=-1.49740e-3,
            d=2.45328e-6,
        ),
    )

    DMC: Component = Component(
        name="DMC",
        molecular_weight=90.08,
        vapour_pressure_constants=VaporPressureConstants(
            a=5.78894,
            b=-1049.375,
            c=-85.977,
        ),
        heat_capacity_constants=HeatCapacityConstants(
            a=3.29155e2,
            b=-1.54438,
            c=4.42614e-3,
            d=-3.66196e-6,
        ),
    )

    CycloHexane: Component = Component(
        name="CycloHexane",
        molecular_weight=84.16,
        vapour_pressure_constants=VaporPressureConstants(
            a=5.97636182,
            b=-1206.47,
            c=-49.864,
        ),
        heat_capacity_constants=HeatCapacityConstants(
            a=43.37101,
            b=0.380403,
            c=0,
            d=0,
        ),
    )

    Benzene: Component = Component(
        name="Benzene",
        molecular_weight=78.11,
        vapour_pressure_constants=VaporPressureConstants(
            a=6.00477182,
            b=-1196.76,
            c=-53.839,
        ),
        heat_capacity_constants=HeatCapacityConstants(
            a=1.18680e2,
            b=-1.01465e-1,
            c=5.41068e-4,
            d=0,
        ),
    )

    Toluene: Component = Component(
        name="Toluene",
        molecular_weight=92.14,
        vapour_pressure_constants=VaporPressureConstants(
            a=6.0854,
            b=-1348.77,
            c=-53.024,
        ),
        heat_capacity_constants=HeatCapacityConstants(
            a=1.87438e2,
            b=-7.30265e-1,
            c=2.96136e-3,
            d=-2.86617e-6,
        ),
    )

    AceticAcid: Component = Component(
        name="Acetic acid",
        molecular_weight=60.05,
        vapour_pressure_constants=VaporPressureConstants(
            a=6.68206,
            b=-1642.54,
            c=-39.764,
        ),
        heat_capacity_constants=HeatCapacityConstants(
            a=-3.891576e3,
            b=3.640888e1,
            c=-1.099113e-1,
            d=1.1060630e-4,
        ),
    )
