from component import AllComponents
from conditions import Conditions
from experiments import IdealExperiments, IdealExperiment
from membrane import Membrane
from mixture import (
    AllMixtures,
    Composition,
    CompositionType,
    get_nrtl_partial_pressures,
)
from pervaporation import Pervaporation

all_components: AllComponents = AllComponents.load("../../components.yml")
all_mixtures: AllMixtures = AllMixtures.load("../../mixtures.yml", all_components)

experiment_meoh_1 = IdealExperiment(
    name="PAI_SPI_1(wt%)_asym",
    temperature=293.15,
    component=all_components.meoh,
    permeance=0.0267,
    activation_energy=20904.65064,
)

experiment_meoh_2 = IdealExperiment(
    name="PAI_SPI_1(wt%)_asym",
    temperature=313.15,
    component=all_components.meoh,
    permeance=0.04498,
    activation_energy=20904.65064,
)

experiment_meoh_3 = IdealExperiment(
    name="PAI_SPI_1(wt%)_asym",
    temperature=325.15,
    component=all_components.meoh,
    permeance=0.06246,
    activation_energy=20904.65064,
)

experiment_mtbe_1 = IdealExperiment(
    name="PAI_SPI_1(wt%)_asym",
    temperature=293.15,
    component=all_components.mtbe,
    permeance=0.01077,
    activation_energy=-39656.76576,
)

experiment_mtbe_2 = IdealExperiment(
    name="PAI_SPI_1(wt%)_asym",
    temperature=313.15,
    component=all_components.mtbe,
    permeance=0.00425,
    activation_energy=-39656.76576,
)

experiment_mtbe_3 = IdealExperiment(
    name="PAI_SPI_1(wt%)_asym",
    temperature=325.15,
    component=all_components.mtbe,
    permeance=0.00212,
    activation_energy=-39656.76576,
)

ideal_experiments = IdealExperiments(
    experiments=[
        experiment_meoh_1,
        experiment_meoh_2,
        # experiment_meoh_3,
        experiment_mtbe_1,
        experiment_mtbe_2,
        # experiment_mtbe_3,
    ]
)

pai_spi = Membrane(ideal_experiments=ideal_experiments, name="PAI_SPI_1(wt%)_asym")

test_conditions = Conditions(
    membrane_area=0.05,
    feed_temperature=323.15,
    permeate_temperature=1,
    feed_amount=1,
    initial_feed_composition=Composition(p=0.15, type=CompositionType("weight")),
)
meoh_mtbe_pervaporation = Pervaporation(
    membrane=pai_spi, mixture=all_mixtures.meoh_mtbe
)

# Check Activation Energy values
activiation_energy_methanol = pai_spi.calculate_activation_energy(all_components.meoh)
activation_energy_mtbe = pai_spi.calculate_activation_energy(all_components.mtbe)
print(
    "Activation Energies: methanol: ",
    activiation_energy_methanol,
    "J/mol  mtbe:",
    activation_energy_mtbe,
    "J/mol",
)

feed_compositions = [
    Composition(p=(0.15 - i / 100), type=CompositionType("weight")) for i in range(15)
]

# Validation of Diffusion curve according to #DOI: 10.1002/app.49982
for i in range(len(feed_compositions)):
    partial_fluxes = meoh_mtbe_pervaporation.calculate_partial_fluxes(
        feed_temperature=325.45, composition=feed_compositions[i]
    )
    psat = get_nrtl_partial_pressures(
        333.15, all_mixtures.meoh_mtbe, feed_compositions[i]
    )
    print(
        partial_fluxes[0],
        " ",
        partial_fluxes[1],
    )

# Check functions of precision
precision_range = [1, 0.5, 0.05, 0.005, 0.0005, 0.00005, 5e-10]

for i in range(len(precision_range)):
    partial_fluxes = meoh_mtbe_pervaporation.calculate_partial_fluxes(
        feed_temperature=333.35,
        composition=Composition(p=0.03, type=CompositionType("weight")),
        precision=precision_range[i],
        permeate_temperature=298,
    )
    print(
        "Precision ",
        precision_range[i],
        " MeOH Flux:",
        partial_fluxes[0],
        " MTBE Flux:",
        partial_fluxes[1],
        " MTBE in permeate:",
        partial_fluxes[1] / sum(partial_fluxes),
    )
