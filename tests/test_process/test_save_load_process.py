import shutil
from pathlib import Path

from pyvaporation.conditions import Conditions
from pyvaporation.membrane import Membrane
from pyvaporation.mixtures import Composition, CompositionType, Mixtures
from pyvaporation.permeance import Permeance
from pyvaporation.pervaporation import Pervaporation
from pyvaporation.process import ProcessModel


def test_save_load_process():

    pervap = Membrane.load(Path("tests/default_membranes/Pervap_4101"))
    romakon = Membrane.load(Path("tests/default_membranes/RomakonPM_102"))

    pv_pervap = Pervaporation(
        membrane=pervap,
        mixture=Mixtures.H2O_EtOH,
    )

    pv_romakon = Pervaporation(
        membrane=romakon,
        mixture=Mixtures.H2O_EtOH,
    )

    con = Conditions(
        membrane_area=0.017,
        initial_feed_temperature=368.15,
        initial_feed_amount=1.5,
        initial_feed_composition=Composition(p=0.1, type=CompositionType.weight),
        permeate_pressure=0,
    )

    non_ideal_process = pv_pervap.non_ideal_isothermal_process(
        conditions=con,
        diffusion_curve_set=pervap.diffusion_curve_sets[0],
        initial_permeances=(
            Permeance(0.0153),
            Permeance(0.00000632),
        ),
        number_of_steps=50,
        delta_hours=0.2,
    )

    ideal_process = pv_romakon.ideal_isothermal_process(
        conditions=con,
        number_of_steps=50,
        delta_hours=0.2,
    )

    process_list = [ideal_process, non_ideal_process]

    non_ideal_process.membrane_path = Path("tests/default_membranes/Pervap_4101")
    ideal_process.membrane_path = Path("tests/default_membranes/RomakonPM_102")

    truth_list = [True, False]

    for process in process_list:
        for is_safe in truth_list:
            process.save(process.membrane_path, is_safe)
            results_path = process.membrane_path / "results"
            process_path = list(
                filter(
                    lambda x: x.stem.startswith("process"),
                    results_path.iterdir(),
                )
            )[0]
            loaded = ProcessModel.load(process_path=process_path, is_safe=is_safe)
            shutil.rmtree(process_path)

            for i in range(len(loaded.time)):
                assert round(loaded.time[i], 2) == round(process.time[i], 2)
                assert round(loaded.partial_fluxes[i][0], 4) == round(
                    process.partial_fluxes[i][0], 4
                )
                assert round(loaded.partial_fluxes[i][0], 4) == round(
                    process.partial_fluxes[i][0], 4
                )
                assert round(loaded.feed_mass[i], 4) == round(process.feed_mass[i], 4)
