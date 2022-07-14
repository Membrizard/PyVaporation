from ..components import Components
from ..utils import NRTLParameters
from .mixture import Mixture


class Mixtures:
    """
    A class with pre-defined mixtures with constants validated against literature data
    For details see ./tests/test_mixtures/test_default_mixtures.py
    """

    H2O_MeOH: Mixture = Mixture(
        name="H2O_MeOH",
        first_component=Components.H2O,
        second_component=Components.MeOH,
        nrtl_params=NRTLParameters(
            g12=-5132.51739,
            g21=1438.40193,
            alpha12=0,
            alpha21=0.3,
            a12=2.7321,
            a21=-0.693,
        ),
    )

    H2O_EtOH: Mixture = Mixture(
        name="H2O_EtOH",
        first_component=Components.H2O,
        second_component=Components.EtOH,
        nrtl_params=NRTLParameters(
            g12=5823,
            g21=-633,
            alpha12=0.3,
        ),
    )

    H2O_iPOH: Mixture = Mixture(
        name="H2O_iPOH",
        first_component=Components.H2O,
        second_component=Components.iPOH,
        nrtl_params=NRTLParameters(
            g12=6899.21,
            g21=106.99,
            alpha12=0.3,
        ),
    )

    H2O_AceticAcid: Mixture = Mixture(
        name="H2O_AceticAcid",
        first_component=Components.H2O,
        second_component=Components.AceticAcid,
        nrtl_params=NRTLParameters(
            g12=-352.42,
            g21=715.43,
            alpha12=0.25,
        ),
    )

    EtOH_ETBE: Mixture = Mixture(
        name="EtOH_ETBE",
        first_component=Components.EtOH,
        second_component=Components.ETBE,
        nrtl_params=NRTLParameters(
            g12=1140.7722,
            g21=2069.17502,
            alpha12=0.3,
        ),
    )

    MeOH_Toluene: Mixture = Mixture(
        name="MeOH_Toluene",
        first_component=Components.MeOH,
        second_component=Components.Toluene,
        nrtl_params=NRTLParameters(
            g12=3857.3,
            g21=4290.3,
            alpha12=0.4370,
        ),
    )

    MeOH_MTBE: Mixture = Mixture(
        name="MeOH_MTBE",
        first_component=Components.MeOH,
        second_component=Components.MTBE,
        nrtl_params=NRTLParameters(
            g12=2133.5,
            g21=2025.3,
            alpha12=0.6,
        ),
    )

    MeOH_DMC: Mixture = Mixture(
        name="MeOH_DMC",
        first_component=Components.MeOH,
        second_component=Components.DMC,
        nrtl_params=NRTLParameters(
            g12=3115.2,
            g21=833.1,
            alpha12=0.3,
        ),
    )
