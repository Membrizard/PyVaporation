from ..components import Components
from ..utils import (
    NRTLParameters,
    UNIQUACParameters,
    UNIQUACBinaryInteractionParameters,
)
from .mixture import Mixture


class Mixtures:
    """
    A class with pre-defined mixtures with constants validated against literature data
    For details see ./tests/test_mixtures/test_default_mixtures.py
    """
    H2O_MeOH: Mixture = Mixture(
        name="H2O_MeOH",
        components=[Components.H2O, Components.MeOH],
        nrtl_params=NRTLParameters(
            g12=-5132.51739,
            g21=1438.40193,
            alpha12=0,
            alpha21=0.3,
            a12=2.7321,
            a21=-0.693,
        ),
        uniquac_params=UNIQUACParameters(
            binary_parameters=[
                UNIQUACBinaryInteractionParameters(
                    i_component_name="H2O",
                    j_component_name="MeOH",
                    ij_parameter=(-3.983130612278737, 0.007057570080263636),
                    ji_parameter=(4.87905141883045, 0.07115788640105822),
                )
            ],
            z=10,
        ),
    )

    H2O_EtOH: Mixture = Mixture(
        name="H2O_EtOH",
        components=[Components.H2O, Components.EtOH],
        nrtl_params=NRTLParameters(
            g12=5823,
            g21=-633,
            alpha12=0.3,
        ),
        uniquac_params=UNIQUACParameters(
            binary_parameters=[
                UNIQUACBinaryInteractionParameters(
                    i_component_name="H2O",
                    j_component_name="EtOH",
                    ij_parameter=(21.127561704493143, -0.9175664931087569),
                    ji_parameter=(100.10268878024358, 2.4619377106475753),
                )
            ],
            z=13,
        ),
    )

    H2O_iPOH: Mixture = Mixture(
        name="H2O_iPOH",
        components=[Components.H2O, Components.iPOH],
        nrtl_params=NRTLParameters(
            g12=6899.21,
            g21=106.99,
            alpha12=0.3,
        ),
        uniquac_params=UNIQUACParameters(
            binary_parameters=[
                UNIQUACBinaryInteractionParameters(
                    i_component_name=Components.H2O.name,
                    j_component_name=Components.iPOH.name,
                    ij_parameter=(-756.0564435691869, 333788.20320719096),
                    ji_parameter=(-678.5581809818217, 266824.57316704467),
                )
            ],
            z=10,
        ),
    )

    H2O_AceticAcid: Mixture = Mixture(
        # TODO: Fix UNIQUAC coefficients to fit with precision in tests
        name="H2O_AceticAcid",
        components=[Components.H2O, Components.AceticAcid],
        nrtl_params=NRTLParameters(
            g12=-352.42,
            g21=715.43,
            alpha12=0.25,
        ),
        uniquac_params=UNIQUACParameters(
            binary_parameters=[
                UNIQUACBinaryInteractionParameters(
                    i_component_name=Components.H2O.name,
                    j_component_name=Components.AceticAcid.name,
                    ij_parameter=(-129.9765770340468, -0.31354951657408175),
                    ji_parameter=(23.640931620853934, 0.07553464090651879),
                )
            ],
            z=10,
        ),
    )

    EtOH_ETBE: Mixture = Mixture(
        # TODO: Fix UNIQUAC coefficients to fit with precision in tests
        name="EtOH_ETBE",
        components=[Components.EtOH, Components.ETBE],
        nrtl_params=NRTLParameters(
            g12=1140.7722,
            g21=2069.17502,
            alpha12=0.3,
        ),
        uniquac_params=UNIQUACParameters(
            binary_parameters=[
                UNIQUACBinaryInteractionParameters(
                    i_component_name=Components.EtOH.name,
                    j_component_name=Components.ETBE.name,
                    ij_parameter=(8942.589565297398, -2915794.8390840776),
                    ji_parameter=(-18441.776691222403, 6478988.764315034),
                )
            ],
            z=10,
        ),
    )

    MeOH_Toluene: Mixture = Mixture(
        name="MeOH_Toluene",
        components=[Components.MeOH, Components.Toluene],
        nrtl_params=NRTLParameters(
            g12=3857.3,
            g21=4290.3,
            alpha12=0.4370,
        ),
        uniquac_params=UNIQUACParameters(
            binary_parameters=[
                UNIQUACBinaryInteractionParameters(
                    i_component_name=Components.MeOH.name,
                    j_component_name=Components.Toluene.name,
                    ij_parameter=(4686.596023943361, -1234765.6427427588),
                    ji_parameter=(-2095.043872895277, 957134.982560884),
                )
            ],
            z=10,
        ),
    )

    MeOH_MTBE: Mixture = Mixture(
        # TODO: Fix UNIQUAC coefficients to fit with precision in tests
        name="MeOH_MTBE",
        components=[Components.MeOH, Components.MTBE],
        nrtl_params=NRTLParameters(
            g12=2133.5,
            g21=2025.3,
            alpha12=0.6,
        ),
        uniquac_params=UNIQUACParameters(
            binary_parameters=[
                UNIQUACBinaryInteractionParameters(
                    i_component_name=Components.MeOH.name,
                    j_component_name=Components.MTBE.name,
                    ij_parameter=(-2487.7680701255767, 903707.2728090351),
                    ji_parameter=(-1614.4771614656215, 651311.6984954888),
                )
            ],
            z=10,
        ),
    )

    MeOH_DMC: Mixture = Mixture(
        # TODO: Fix UNIQUAC coefficients to fit with precision in tests
        name="MeOH_DMC",
        components=[Components.MeOH, Components.DMC],
        nrtl_params=NRTLParameters(
            g12=3115.2,
            g21=833.1,
            alpha12=0.3,
        ),
        uniquac_params=UNIQUACParameters(
            binary_parameters=[
                UNIQUACBinaryInteractionParameters(
                    i_component_name=Components.MeOH.name,
                    j_component_name=Components.DMC.name,
                    ij_parameter=(739.9268135127102, -173417.54480148194),
                    ji_parameter=(-168.38470470351714, 72635.51155280948),
                )
            ],
            z=10,
        ),
    )
