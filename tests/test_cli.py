from typer.testing import CliRunner

from lanctools.cli import app

runner = CliRunner()


def test_version():
    result = runner.invoke(app, ["version"])

    assert result.exit_code == 0
    assert result.stdout.startswith("lanctools ")


def test_merge_with_repeated_inputs(monkeypatch):
    calls = []

    def fake_merge(files, outfile):
        calls.append((files, outfile))

    monkeypatch.setattr("lanctools.merge_lanc", fake_merge)

    result = runner.invoke(
        app,
        [
            "merge",
            "--input",
            "chr1.lanc",
            "--input",
            "chr2.lanc",
            "--output",
            "merged.lanc",
        ],
    )

    assert result.exit_code == 0
    assert calls == [(["chr1.lanc", "chr2.lanc"], "merged.lanc")]


def test_merge_with_input_list(tmp_path, monkeypatch):
    input_list = tmp_path / "inputs.txt"
    input_list.write_text("chr1.lanc\n\nchr2.lanc\n")
    calls = []

    def fake_merge(files, outfile):
        calls.append((files, outfile))

    monkeypatch.setattr("lanctools.merge_lanc", fake_merge)

    result = runner.invoke(
        app,
        ["merge", "--input-list", str(input_list), "--output", "merged.lanc"],
    )

    assert result.exit_code == 0
    assert calls == [(["chr1.lanc", "chr2.lanc"], "merged.lanc")]


def test_merge_rejects_multiple_input_sources(tmp_path):
    input_list = tmp_path / "inputs.txt"
    input_list.write_text("chr1.lanc\n")

    result = runner.invoke(
        app,
        [
            "merge",
            "--input",
            "chr1.lanc",
            "--input-list",
            str(input_list),
            "--output",
            "merged.lanc",
        ],
    )

    assert result.exit_code == 2
    assert "either input OR input_list" in result.output


def test_merge_requires_input():
    result = runner.invoke(app, ["merge", "--output", "merged.lanc"])

    assert result.exit_code == 2
    assert "Specify one of either input or input_list" in result.output


def test_convert_dispatches_to_core(monkeypatch):
    calls = []

    def fake_convert(**kwargs):
        calls.append(kwargs)

    monkeypatch.setattr("lanctools.convert_to_lanc", fake_convert)

    result = runner.invoke(
        app,
        [
            "convert",
            "--input",
            "chr1.anc.vcf.gz",
            "--plink",
            "chr1",
            "--format",
            "FLARE",
            "--output",
            "chr1.lanc",
        ],
    )

    assert result.exit_code == 0
    assert calls == [
        {
            "file": "chr1.anc.vcf.gz",
            "file_fmt": "FLARE",
            "plink_prefix": "chr1",
            "output": "chr1.lanc",
        }
    ]


def test_convert_propagates_core_errors(monkeypatch):
    def fail_convert(**kwargs):
        raise ValueError("unsupported input")

    monkeypatch.setattr("lanctools.convert_to_lanc", fail_convert)

    result = runner.invoke(
        app,
        [
            "convert",
            "--input",
            "input.txt",
            "--plink",
            "chr1",
            "--format",
            "RFMix",
            "--output",
            "chr1.lanc",
        ],
    )

    assert result.exit_code == 1
    assert isinstance(result.exception, ValueError)
    assert str(result.exception) == "unsupported input"
