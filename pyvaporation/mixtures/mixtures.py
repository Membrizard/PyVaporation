from ..components import Components
from ..utils import NRTLParameters, UNIQUACParameters
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
        uniquac_params=UNIQUACParameters(
            alpha_12=-1.5236194657557556,
            alpha_21=2.4511555703523626,
            beta_12=-0.009196570921123198,
            beta_21=0.041136887744315336,
            z=10,
        )
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

        uniquac_params=UNIQUACParameters(
            alpha_12=21.127561704493143,
            alpha_21=100.10268878024358,
            beta_12=-0.9175664931087569,
            beta_21=2.4619377106475753,
            z=13,
        )
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
        uniquac_params=UNIQUACParameters(
            alpha_12=104.55510898284338,
            alpha_21=8.60861483037599,
            beta_12=-6.7637165811854265,
            beta_21=1.8578246644616325,
            z=16,
        )
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
        uniquac_params=UNIQUACParameters(
            alpha_12=-11.138624089082999,
            alpha_21=3.9488095210545184,
            beta_12=1.9362990236330024,
            beta_21=1.7234361729360623,
            z=6,
        )
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
        uniquac_params=UNIQUACParameters(
            alpha_12=-169.03517903712503,
            alpha_21=-74.56587147924621,
            beta_12=12327.525553727144,
            beta_21=-16498.35261855843,
            z=1064,
        )
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
        uniquac_params=UNIQUACParameters(
            alpha_12=13.252,
            alpha_21=-106.210,
            beta_12=368.586,
            beta_21=-57.820,
            z=10,
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
