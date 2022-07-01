from pathlib import Path

from pytest import fixture

from membrane import Membrane
from mixtures import Composition
from permeance import Permeance


@fixture
def diffusion_curve():

    membrane = Membrane.load(Path("tests/default_membranes/Pervap_2510"))
    return membrane.diffusion_curve_sets[0][0]


def test_permeate_composition(diffusion_curve):
    assert diffusion_curve.permeate_composition == [
        Composition(p=0.9962264150943396, type="weight"),
        Composition(p=0.9949238578680203, type="weight"),
        Composition(p=0.9940209267563528, type="weight"),
        Composition(p=0.9964664310954063, type="weight"),
    ]


def test_get_separation_factor(diffusion_curve):
    validation_separation_factor = [
        2122.5414841475217,
        1298.2177925790918,
        1069.644713092987,
        1607.6269420984008,
    ]

    for i in range(len(diffusion_curve.get_separation_factor)):
        assert abs(diffusion_curve.get_separation_factor[i]-validation_separation_factor[i]) < 1e-2



def test_get_psi(diffusion_curve):
    validation_psi = [1124.4169866, 766.65571541, 714.92331306, 1364.02627384]
    tested_psi = diffusion_curve.get_psi
    for i in range(len(validation_psi)):
        assert round(validation_psi[i], 1) == round(tested_psi[i], 1)


def test_get_permeances(diffusion_curve):

    validation_permeances = [
        (
            Permeance(value=0.038884957564070496, units="kg/(m2*h*kPa)"),
            Permeance(value=6.839607107275804e-05, units="kg/(m2*h*kPa)"),
        ),
        (
            Permeance(value=0.04014278123404922, units="kg/(m2*h*kPa)"),
            Permeance(value=0.00010620585974979516, units="kg/(m2*h*kPa)"),
        ),
        (
            Permeance(value=0.04493732503525243, units="kg/(m2*h*kPa)"),
            Permeance(value=0.00014235012719086482, units="kg/(m2*h*kPa)"),
        ),
        (
            Permeance(value=0.0549547245885328, units="kg/(m2*h*kPa)"),
            Permeance(value=0.00010910826951395723, units="kg/(m2*h*kPa)"),
        ),
    ]

    for i in range(len(diffusion_curve.get_permeances)):
        assert abs(diffusion_curve.get_permeances[i][0].value - validation_permeances[i][0].value) < 1e-6
        assert abs(diffusion_curve.get_permeances[i][1].value - validation_permeances[i][1].value) < 1e-6


def test_get_selectivity(diffusion_curve):
    validation_selectivity = [
        1896.138951240108,
        1260.60395417634,
        1052.8561604828328,
        1679.8366537203292,
    ]

    for i in range(len(diffusion_curve.get_selectivity)):
        assert (diffusion_curve.get_selectivity[i] - validation_selectivity[i]) < 1e-2
