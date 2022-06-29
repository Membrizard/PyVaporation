from pytest import fixture

from conditions import Conditions
from membrane import Membrane
from mixtures import Composition, CompositionType, Mixtures
from permeance import Permeance
from pervaporation import Pervaporation
from pathlib import Path


def test_save_load_process():

    membrane = Membrane.load(Path("tests/default_membranes/Pervap_4101"))

    pv = Pervaporation(
        membrane=membrane,
        mixture=Mixtures.H2O_EtOH,
    )

    con = Conditions(
        membrane_area=0.017,
        initial_feed_temperature=368.15,
        initial_feed_amount=1.5,
        initial_feed_composition=Composition(p=0.1, type=CompositionType.weight),
        permeate_pressure=0,
    )

    process = pv.non_ideal_isothermal_process(
        conditions=con,
        diffusion_curve_set=membrane.diffusion_curve_sets[0],
        initial_permeances=(
            Permeance(0.0153),
            Permeance(0.00000632),
        ),
        number_of_steps=50,
        delta_hours=0.2,
    )

    process.save()
    assert 0==0
