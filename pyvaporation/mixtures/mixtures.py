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
            alpha_12=-3.983130612278737,
            alpha_21=4.87905141883045,
            beta_12=0.007057570080263636,
            beta_21=0.07115788640105822,
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
            alpha_12=-756.0564435691869,
            alpha_21=-678.5581809818217,
            beta_12=333788.20320719096,
            beta_21=266824.57316704467,
            z=10,
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
            alpha_12=-129.9765770340468,
            alpha_21=23.640931620853934,
            beta_12=-0.31354951657408175,
            beta_21=0.07553464090651879,
            z=10,
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
            alpha_12=8942.589565297398,
            alpha_21=-18441.776691222403,
            beta_12=-2915794.8390840776,
            beta_21=6478988.764315034,
            z=10,
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
            alpha_12=4686.596023943361,
            alpha_21=-2095.043872895277,
            beta_12=-1234765.6427427588,
            beta_21=957134.982560884,
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
        uniquac_params=UNIQUACParameters(
            alpha_12=-2487.7680701255767,
            alpha_21=-1614.4771614656215,
            beta_12=903707.2728090351,
            beta_21=651311.6984954888,
            z=10,
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
        uniquac_params=UNIQUACParameters(
            alpha_12=739.9268135127102,
            alpha_21=-168.38470470351714,
            beta_12=-173417.54480148194,
            beta_21=72635.51155280948,
            z=473,
        ),
    )
