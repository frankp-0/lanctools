# MIT License
# Copyright (c) 2025 Franklin Ockerman
# See LICENSE file for full license text

"""The command line interface for lanctools."""

import logging
import typer
from typing import Optional
from importlib.metadata import version, PackageNotFoundError

app = typer.Typer(help="lanctools CLI")
logger = logging.getLogger("lanctools")


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
        return version("lanctools")
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
    typer.echo(f"lanctools {get_version()}")


@app.command()
def merge(
    input: Optional[list[str]] = typer.Option(
        None,
        help=(
            "Local ancestry files to be merged. "
            "This option should be repeated to specify multiple files. "
            "Example: --input chr1.lanc --input chr2.lanc"
        ),
    ),
    input_list: Optional[str] = typer.Option(
        None,
        help=(
            "File containing paths to local ancestry files to be merged, one per line"
        ),
    ),
    output: str = typer.Option(..., help="Path of output Local ancestry file"),
):
    from . import merge_lanc

    if input and input_list:
        raise typer.BadParameter("Specify either input OR input_list, not both.")

    lancs: list[str]
    if input is None:
        if input_list is None:
            raise typer.BadParameter("Specify one of either input or input_list.")
        with open(input_list) as f:
            lancs = [line.strip() for line in f if line.strip()]
    else:
        lancs = input

    merge_lanc(lancs, output)


@app.command()
def convert(
    input: str = typer.Option(..., help="Local ancestry file"),
    plink: str = typer.Option(..., help="Plink2 file prefix"),
    format: str = typer.Option(..., help="File format (either 'RFMix' or 'FLARE'"),
    output: str = typer.Option(
        ...,
        help="Output prefix",
    ),
):
    from . import convert_to_lanc

    convert_to_lanc(file=input, file_fmt=format, plink_prefix=plink, output=output)


def main_entry():
    try:
        app()
    except Exception as exc:
        logger.debug("Unhandled exception", exc_info=True)
        typer.secho(f"Error: {exc}", fg=typer.colors.RED)
        raise typer.Exit(code=1)


if __name__ == "__main__":
    main_entry()
