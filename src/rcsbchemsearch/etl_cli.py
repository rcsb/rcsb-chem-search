import enum
import re
from enum import StrEnum

from dataclasses import dataclass
from functools import cached_property
from typing import Annotated, Literal, ClassVar

from typer import Typer, Argument, Option
from pathlib import Path
from rcsbchemsearch.etl import Pipeline
from loguru import logger

from rcsbchemsearch.core.rdkit_utils import RdkitLogging


# Don't allow TRACE, SUCCESS, and CRITICAL for broader compatibility.
class LogLevel(StrEnum):
    DEBUG = enum.auto()
    INFO = enum.auto()
    WARNING = enum.auto()
    ERROR = enum.auto()


app = Typer()


@dataclass(frozen=True, slots=True)
class LogSink:
    _PATTERN: ClassVar[re.Pattern] = re.compile(r"(?P<pipe>STDOUT|STDERR)|(?P<file>.+?\.log(?:\.gz)?)")
    name: str

    def __post_init__(self) -> None:
        if not self._PATTERN.fullmatch(self.name):
        if self.name not in {"STDOUT", "STDERR"} and not self.name.endswith(".log") and not self.name.endswith(".gz"):
            msg = f"Invalid log sink '{self.name}'."
            raise ValueError(msg)

    @cached_property
    def path(self) -> Path | None:
        return None if self.name in {"STDOUT", "STDERR"} else Path(self.name)

    @cached_property
    def compression(self) -> Literal["gzip"] | None:
        return "gzip" if isinstance(self.path, Path) and self.path.suffix == ".gz" else None


@app.command()
def run(
    input_file: Annotated[
        str,
        Argument(help="URI to the input JSON file containing molecules."),
    ],
    output_path: Annotated[
        Path,
        Argument(help="Path for the output `.json.zst` file."),
    ],
    fingerprints: Annotated[
        frozenset[str],
        Option(help="List of fingerprint types to calculate (comma-separated).", parser=str.split),
    ] = frozenset(),
    max_tautomers: Annotated[
        int,
        Option(help="Maximum number of tautomers to generate per molecule.", min=1),
    ] = 1,
    min_tautomer_score: Annotated[
        int,
        Option(help="Minimum quality score for accepting tautomers.", min=0),
    ] = 0,
    max_tautomer_transforms: Annotated[
        int,
        Option(help="Maximum number of total transforms per molecule to use when generating tautomers.", min=0),
    ] = 1000,
    random_seed: Annotated[
        int,
        Option(help="Random seed for reproducibility.", min=1),
    ] = 1,
    resume: Annotated[
        bool,
        Option(help="Flag to resume processing from the last checkpoint."),
    ] = False,
    halt_on_error: Annotated[
        bool,
        Option(help=""),
    ] = False,
    log_level: Annotated[
        LogLevel,
        Option(help="Minimum log level."),
    ] = LogLevel.INFO,
    log_sink: Annotated[
        LogSink,
        Option(help="STDERR, STDOUT, a file:// URI, or a path ending in .log or .log.gz."),
    ] = str,
):
    """
    Run the ETL pipeline to process molecules and generate fingerprints.
    """
    logger.configure(
        handlers=[]
    )
    Log.apply_config(LogConfigFactory().of(log_level, log_sink))
    Pipeline(
        input_file=input_file,
        output_path=output_path,
        fingerprints=fingerprints,
        max_tautomers=max_tautomers,
        min_tautomer_score=min_tautomer_score,
        resume=resume,
        random_sed=random_sed,
        halt_on_error=halt_on_error,
    ).run()


if __name__ == "__main__":
    app()
