import typing
from pathlib import Path

import attr
import yaml


def _str_to_path(value: str) -> Path:
    return Path(value)


@attr.s(auto_attribs=True)
class Config:
    source_path: Path = attr.ib(converter=_str_to_path)
    # results_path: Path = attr.ib(converter=_str_to_path)

    def __attrs_post_init__(self):
        self.ideal_experiments_path = self.source_path / "ideal_experiments.csv"
        self.diffusion_curve_sets_path = self.source_path / "diffusion_curve_sets"

        if not (
            self.ideal_experiments_path.exists()
            or self.diffusion_curve_sets_path.exists()
        ):
            raise FileExistsError("Missing experiment data at %s" % self.source_path)

    @classmethod
    def load(cls, path: typing.Union[str, Path]) -> "Config":
        with open(path, "rb") as handle:
            return cls(**yaml.load(handle, Loader=yaml.FullLoader))
