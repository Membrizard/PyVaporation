from component import AllComponents
from conditions import Conditions
from experiments import IdealExperiment, IdealExperiments
from membrane import Membrane
from mixture import (AllMixtures, Composition, CompositionType,
                     get_nrtl_partial_pressures)
from pervaporation import Pervaporation

all_components: AllComponents = AllComponents.load("../../components.yml")
all_mixtures: AllMixtures = AllMixtures.load("../../mixtures.yml", all_components)

experiment_h2o_1 = IdealExperiment(
    name="Pervap 4101",
    temperature=333.15,
    component=all_components.h2o,
    permeance=0.01491,
    activation_energy=65156.28,
)

experiment_h2o_2 = IdealExperiment(
    name="Pervap 4101",
    temperature=368.15,
    component=all_components.h2o,
    permeance=0.02460,
    activation_energy=65156.28,
)

experiment_h2o_3 = IdealExperiment(
    name="Pervap 4101",
    temperature=378.15,
    component=all_components.h2o,
    permeance=0.03555,
    activation_energy=65156.28,
)

experiment_etoh_1 = IdealExperiment(
    name="Pervap 4101",
    temperature=368.15,
    component=all_components.etoh,
    permeance=0.0000098,
)

ideal_experiments = IdealExperiments(
    experiments=[
        experiment_h2o_1,
        experiment_h2o_2,
        experiment_h2o_3,
        experiment_etoh_1,
    ]
)

pervap_4101 = Membrane(ideal_experiments=ideal_experiments, name="Pervap 4101")

test_conditions = Conditions(
    membrane_area=0.05,
    feed_temperature=323.15,
    permeate_temperature=1,
    feed_amount=1,
    initial_feed_composition=Composition(p=0.15, type=CompositionType("weight")),
)
h2o_etoh_pervaporation = Pervaporation(
    membrane=pervap_4101, mixture=all_mixtures.h2o_etoh
)

# Check Activation Energy values
activiation_energy_h2o = pervap_4101.calculate_activation_energy(all_components.h2o)
activation_energy_etoh = 0
print(
    "Activation Energies: h2o: ",
    activiation_energy_h2o,
    "J/mol  etoh:",
    activation_energy_etoh,
    "J/mol",
)

feed = [
    0.0080,
    0.0105,
    0.0193,
    0.0247,
    0.0315,
    0.0585,
    0.0728,
    0.0869,
    0.1000,
    0.1341,
    0.1536,
]

feed_compositions = [
    Composition(p=feed[i], type=CompositionType("weight")) for i in range(len(feed))
]

# Validation of Diffusion curve according to DOI: 10.22079/JMSR.2018.88186.1198
for i in range(len(feed_compositions)):
    partial_fluxes = h2o_etoh_pervaporation.calculate_partial_fluxes(
        feed_temperature=368.15, composition=feed_compositions[i]
    )
    psat = get_nrtl_partial_pressures(
        368.15, all_mixtures.h2o_etoh, feed_compositions[i]
    )
    print(
        "EtOH in feed",
        feed_compositions[i].second,
        " Total Flux ",
        sum(partial_fluxes),
        " EtOH in Permeate ",
        (partial_fluxes[1] / sum(partial_fluxes)),
        "Psat H2O ",
        psat[0],
    )
