# MIT License
# Copyright (c) 2025 Franklin Ockerman
# See LICENSE file for full license text

import logging
import typer
from typing import Optional, List
from importlib.metadata import version, PackageNotFoundError

app = typer.Typer(help="lagga CLI")
logger = logging.getLogger("lagga")


def list_from_csv(arg: str) -> List[str]:
    return [x.strip() for x in arg.split(",")]


def setup_logging(verbose: bool, quiet: bool) -> None:
    if quiet:
        level = logging.ERROR
    elif verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def get_version() -> str:
    try:
        return version("lanctools")  # must match your project name
    except PackageNotFoundError:
        return "unknown"


@app.callback()
def main(
    verbose: bool = typer.Option(False, "--verbose"),
    quiet: bool = typer.Option(False, "--quiet"),
):
    setup_logging(verbose, quiet)


@app.command(name="version", help="Show the CLI version and exit")
def show_version():
    typer.echo(f"lagga {get_version()}")


@app.command()
def convert_flare(
    file: str = typer.Option(..., help="Local ancestry file(s), comma-separated"),
    plink_prefix: str = typer.Option(
        ..., help="Plink2 file prefix(es), comma-separated"
    ),
    output: str = typer.Option(
        ...,
        help="Output prefix(es), comma-separated, one per plink_prefix",
    ),
):
    from . import convert_to_lanc

    plinks = list_from_csv(plink_prefix)
    inputs = list_from_csv(file)
    outputs = list_from_csv(output)

    for plink, input, out in zip(plinks, inputs, outputs):
        convert_to_lanc(file=input, file_fmt="FLARE", plink_prefix=plink, output=out)


@app.command()
def convert_rfmix(
    file: str = typer.Option(..., help="Local ancestry file(s), comma-separated"),
    plink_prefix: str = typer.Option(
        ..., help="Plink2 file prefix(es), comma-separated"
    ),
    outputs: str = typer.Option(
        ...,
        help="Output prefix(es), comma-separated, one per plink_prefix",
    ),
):
    from . import convert_to_lanc

    plinks = list_from_csv(plink_prefix)
    inputs = list_from_csv(file)
    output = list_from_csv(outputs)

    for plink, input, out in zip(plinks, inputs, output):
        convert_to_lanc(file=input, file_fmt="RFMix", plink_prefix=plink, output=out)


def main_entry():
    try:
        app()
    except Exception as exc:
        logger.debug("Unhandled exception", exc_info=True)
        typer.secho(f"Error: {exc}", fg=typer.colors.RED)
        raise typer.Exit(code=1)


if __name__ == "__main__":
    main_entry()
