import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from lanctools import FlatLanc, LancData, convert_to_lanc, core, merge_lanc
from lanctools._cpp import read_flare, read_rfmix
from lanctools.core import _parse_lanc_line, _read_lanc


@pytest.fixture
def chr20_data():
    dataset = LancData(plink_prefix="tests/data/chr20", lanc_file="tests/data/chr20.lanc")
    return dataset


def test_parse_lanc_basic():
    line = "10:02 14:12 30:11"
    left, right, bp = _parse_lanc_line(line)

    np.testing.assert_array_equal(left, [0, 1, 1])
    np.testing.assert_array_equal(right, [2, 2, 1])
    np.testing.assert_array_equal(bp, [10, 14, 30])


def test_parse_lanc_hom():
    line = "10:00"
    left, right, bp = _parse_lanc_line(line)

    np.testing.assert_array_equal(left, [0])
    np.testing.assert_array_equal(right, [0])
    np.testing.assert_array_equal(bp, [10])


def test_read_lanc(tmp_path):
    content = "6 2\n10:01 20:12 30:23\n15:10 40:01\n"
    path = tmp_path / "test.lanc"
    path.write_text(content)

    lanc = _read_lanc(path)
    np.testing.assert_array_equal(lanc.offsets, [0, 3, 5])
    np.testing.assert_array_equal(lanc.breakpoints, [10, 20, 30, 15, 40])
    np.testing.assert_array_equal(lanc.left_haps, [0, 1, 2, 1, 0])
    np.testing.assert_array_equal(lanc.right_haps, [1, 2, 3, 0, 1])


def test_parse_lanc_rejects_invalid_tract():
    with pytest.raises(ValueError, match="Breakpoints"):
        _parse_lanc_line("10:01 9:12")


def test_read_lanc_rejects_wrong_sample_count(tmp_path):
    path = tmp_path / "invalid.lanc"
    path.write_text("6 2\n6:00\n")

    with pytest.raises(ValueError, match="declares 2 samples"):
        _read_lanc(path)


def test_merge_lanc_rejects_empty_input(tmp_path):
    with pytest.raises(ValueError, match="At least one"):
        merge_lanc([], tmp_path / "output.lanc")


def test_merge_lanc_writes_offset_breakpoints(tmp_path):
    first = tmp_path / "first.lanc"
    first.write_text("3 2\n1:00 3:11\n3:01\n")
    second = tmp_path / "second.lanc"
    second.write_text("2 2\n2:10\n2:11\n")
    output = tmp_path / "merged.lanc"

    merge_lanc([str(first), str(second)], output)

    assert output.read_text() == "5 2\n1:00 3:11 5:10\n3:01 5:11\n"


def test_merge_lanc_rejects_incorrect_final_breakpoint(tmp_path):
    source = tmp_path / "invalid.lanc"
    source.write_text("3 1\n2:00\n")

    with pytest.raises(ValueError, match="final breakpoint"):
        merge_lanc([str(source)], tmp_path / "output.lanc")


def test_merge_lanc_rejects_breakpoint_collision(tmp_path):
    first = tmp_path / "first.lanc"
    first.write_text("3 1\n3:00\n")
    second = tmp_path / "second.lanc"
    second.write_text("1 1\n0:11 1:00\n")

    with pytest.raises(ValueError, match="strictly increasing"):
        merge_lanc([str(first), str(second)], tmp_path / "output.lanc")


def test_merge_lanc_rejects_extra_sample_lines(tmp_path):
    source = tmp_path / "invalid.lanc"
    source.write_text("1 1\n1:00\n1:00\n")

    with pytest.raises(ValueError, match="more sample lines"):
        merge_lanc([str(source)], tmp_path / "output.lanc")


def test_convert_flare(tmp_path):
    tmp_lanc_path = tmp_path / "test_flare.lanc"
    convert_to_lanc(
        file="tests/data/chr20.flare.anc.vcf.gz",
        file_fmt="FLARE",
        plink_prefix="tests/data/chr20",
        output=tmp_lanc_path,
    )

    with (
        open("tests/data/chr20.lanc", encoding="utf-8") as true_lanc,
        open(tmp_lanc_path, encoding="utf-8") as test_lanc,
    ):
        assert true_lanc.read() == test_lanc.read()


def test_convert_rejects_missing_parser_columns(monkeypatch):
    monkeypatch.setattr(core, "read_flare", lambda _: {})

    with pytest.raises(ValueError, match="missing columns"):
        convert_to_lanc(
            file="input.vcf.gz",
            file_fmt="FLARE",
            plink_prefix="tests/data/chr20",
            output="output.lanc",
        )


def test_convert_rejects_empty_parser_output(monkeypatch):
    columns = ["sample", "chrom", "spos", "epos", "anc0", "anc1"]
    monkeypatch.setattr(core, "read_flare", lambda _: {column: [] for column in columns})

    with pytest.raises(ValueError, match="no ancestry tracts"):
        convert_to_lanc(
            file="input.vcf.gz",
            file_fmt="FLARE",
            plink_prefix="tests/data/chr20",
            output="output.lanc",
        )


def test_convert_rejects_unknown_chromosome(monkeypatch):
    original_read_flare = core.read_flare
    flare_data = original_read_flare("tests/data/chr20.flare.anc.vcf.gz")
    flare_data["chrom"] = ["unknown"] * len(flare_data["chrom"])
    monkeypatch.setattr(core, "read_flare", lambda _: flare_data)

    with pytest.raises(ValueError, match="unknown chromosomes"):
        convert_to_lanc(
            file="input.vcf.gz",
            file_fmt="FLARE",
            plink_prefix="tests/data/chr20",
            output="output.lanc",
        )


def test_convert_rejects_sample_mismatch(monkeypatch):
    original_read_flare = core.read_flare
    flare_data = original_read_flare("tests/data/chr20.flare.anc.vcf.gz")
    flare_data["sample"] = ["unexpected"] * len(flare_data["sample"])
    monkeypatch.setattr(core, "read_flare", lambda _: flare_data)

    with pytest.raises(ValueError, match="samples differ"):
        convert_to_lanc(
            file="input.vcf.gz",
            file_fmt="FLARE",
            plink_prefix="tests/data/chr20",
            output="output.lanc",
        )


def test_read_flare_rejects_unreadable_file(tmp_path):
    with pytest.raises(RuntimeError, match="Failed to open input VCF"):
        read_flare(str(tmp_path / "missing.vcf.gz"))


def test_read_flare_rejects_missing_header(tmp_path):
    path = tmp_path / "missing_header.vcf.gz"
    with gzip.open(path, "wt") as vcf:
        vcf.write("##fileformat=VCFv4.3\n")

    with pytest.raises(RuntimeError, match="Missing #CHROM header"):
        read_flare(str(path))


def test_read_rfmix_rejects_odd_haplotype_count(tmp_path):
    path = tmp_path / "odd_haplotypes.msp.tsv"
    path.write_text("population codes\nchrom\tstart\tend\tcm\tn_snps\t...\tsample.0\n")

    with pytest.raises(RuntimeError, match="must be even"):
        read_rfmix(str(path))


def test_read_rfmix_rejects_malformed_record(tmp_path):
    path = tmp_path / "malformed.msp.tsv"
    path.write_text(
        "population codes\nchrom\tstart\tend\tcm\tn_snps\t...\tsample.0\tsample.1\nchr1\t0\t10\n"
    )

    with pytest.raises(RuntimeError, match="too few fields"):
        read_rfmix(str(path))


def test_get_info(chr20_data):
    nvar = chr20_data.pvar.get_variant_ct()
    df_info = chr20_data.get_info(np.arange(nvar))
    df_true = pd.read_json(
        "tests/data/chr20_info.json",
        dtype={"CHR": str, "BP": np.uint32, "REF": str, "ALT": str, "ID": str},
    )
    pd.testing.assert_frame_equal(df_info, df_true)


def test_get_lanc(chr20_data):
    lanc_arr = chr20_data.get_lanc(np.arange(10, 14, dtype=np.uint32))
    lanc_true = np.asarray(
        [
            [[1, 0], [1, 0], [1, 0], [1, 0]],
            [[1, 1], [1, 1], [1, 1], [1, 1]],
            [[1, 1], [1, 1], [1, 1], [1, 1]],
            [[1, 1], [1, 1], [1, 1], [1, 1]],
            [[1, 1], [1, 1], [1, 1], [1, 1]],
            [[0, 1], [0, 1], [0, 1], [0, 1]],
            [[1, 1], [1, 1], [1, 1], [1, 1]],
            [[0, 1], [0, 1], [0, 1], [0, 1]],
            [[0, 1], [0, 1], [0, 1], [1, 1]],
            [[1, 0], [1, 0], [1, 0], [1, 0]],
            [[1, 0], [1, 0], [1, 0], [1, 0]],
            [[1, 1], [1, 1], [1, 1], [1, 1]],
            [[1, 0], [1, 0], [1, 0], [1, 1]],
            [[1, 0], [1, 0], [1, 0], [1, 0]],
            [[1, 0], [1, 0], [1, 0], [1, 0]],
            [[0, 1], [0, 1], [0, 1], [0, 1]],
            [[0, 0], [0, 0], [0, 1], [0, 1]],
            [[1, 1], [1, 0], [1, 0], [1, 0]],
            [[1, 1], [1, 1], [1, 1], [1, 1]],
            [[0, 1], [0, 1], [0, 1], [0, 1]],
        ],
        dtype=np.uint8,
    )

    np.testing.assert_equal(lanc_arr, lanc_true)


def test_get_lanc_rejects_out_of_bounds_indices(chr20_data):
    with pytest.raises(IndexError, match="outside"):
        chr20_data.get_lanc(np.array([40], dtype=np.uint32))


@pytest.mark.parametrize(
    ("indices", "error", "message"),
    [
        (np.array([[1]], dtype=np.uint32), ValueError, "one-dimensional"),
        (np.array([1.0]), TypeError, "integer dtype"),
        (np.array([-1], dtype=np.int32), IndexError, "non-negative"),
        (np.array([40], dtype=np.uint32), IndexError, "outside"),
    ],
)
def test_query_methods_validate_indices(chr20_data, indices, error, message):
    with pytest.raises(error, match=message):
        chr20_data.get_info(indices)

    with pytest.raises(error, match=message):
        chr20_data.get_geno(indices)


def test_lanc_data_rejects_variant_count_mismatch(tmp_path):
    lanc_path = tmp_path / "wrong_variant_count.lanc"
    lines = Path("tests/data/chr20.lanc").read_text().splitlines()
    lanc_path.write_text("39 20\n" + "\n".join(lines[1:]) + "\n")

    with pytest.raises(ValueError, match="numbers of variants"):
        LancData(plink_prefix="tests/data/chr20", lanc_file=str(lanc_path))


def test_lanc_data_rejects_incomplete_ancestry_names():
    with pytest.raises(ValueError, match="Ancestry names"):
        LancData(
            plink_prefix="tests/data/chr20",
            lanc_file="tests/data/chr20.lanc",
            ancestries=["only"],
        )


def test_flat_lanc_rejects_mismatched_arrays():
    with pytest.raises(ValueError, match="equal lengths"):
        FlatLanc(
            np.array([0], dtype=np.uint8),
            np.array([], dtype=np.uint8),
            np.array([1], dtype=np.uint32),
            np.array([0, 1], dtype=np.uint32),
        )


def test_get_lanc_dosage(chr20_data):
    lanc_arr = chr20_data.get_lanc_dosage(np.arange(10, 14, dtype=np.uint32))
    lanc_true = np.asarray(
        [
            [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0]],
            [[0.0, 2.0], [0.0, 2.0], [0.0, 2.0], [0.0, 2.0]],
            [[0.0, 2.0], [0.0, 2.0], [0.0, 2.0], [0.0, 2.0]],
            [[0.0, 2.0], [0.0, 2.0], [0.0, 2.0], [0.0, 2.0]],
            [[0.0, 2.0], [0.0, 2.0], [0.0, 2.0], [0.0, 2.0]],
            [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0]],
            [[0.0, 2.0], [0.0, 2.0], [0.0, 2.0], [0.0, 2.0]],
            [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0]],
            [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [0.0, 2.0]],
            [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0]],
            [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0]],
            [[0.0, 2.0], [0.0, 2.0], [0.0, 2.0], [0.0, 2.0]],
            [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [0.0, 2.0]],
            [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0]],
            [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0]],
            [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0]],
            [[2.0, 0.0], [2.0, 0.0], [1.0, 1.0], [1.0, 1.0]],
            [[0.0, 2.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0]],
            [[0.0, 2.0], [0.0, 2.0], [0.0, 2.0], [0.0, 2.0]],
            [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0]],
        ],
        dtype=np.int32,
    )

    np.testing.assert_equal(lanc_arr, lanc_true)


def test_get_geno(chr20_data):
    geno_arr = chr20_data.get_geno(np.arange(10, 14, dtype=np.uint32))
    geno_true = np.asarray(
        [
            [[0, 0], [0, 0], [0, 0], [0, 0]],
            [[0, 0], [0, 0], [0, 1], [0, 0]],
            [[0, 0], [0, 0], [1, 0], [0, 0]],
            [[0, 0], [0, 0], [0, 1], [0, 0]],
            [[0, 0], [0, 1], [0, 0], [0, 0]],
            [[0, 0], [1, 0], [0, 0], [0, 0]],
            [[0, 1], [0, 0], [0, 0], [0, 0]],
            [[0, 0], [0, 0], [1, 1], [1, 0]],
            [[0, 0], [0, 0], [1, 1], [0, 0]],
            [[0, 0], [0, 1], [0, 1], [0, 0]],
            [[0, 0], [0, 0], [0, 0], [0, 0]],
            [[0, 1], [0, 0], [1, 0], [0, 0]],
            [[0, 0], [0, 0], [0, 0], [0, 0]],
            [[0, 0], [0, 1], [0, 0], [0, 0]],
            [[0, 0], [0, 0], [1, 0], [0, 0]],
            [[0, 0], [0, 0], [0, 0], [1, 0]],
            [[0, 0], [0, 1], [1, 0], [0, 0]],
            [[0, 0], [0, 0], [0, 0], [0, 0]],
            [[0, 0], [1, 0], [1, 1], [0, 0]],
            [[0, 0], [0, 0], [0, 0], [0, 0]],
        ],
        dtype=np.int32,
    )

    np.testing.assert_equal(geno_arr, geno_true)


def test_get_lanc_geno(chr20_data):
    lanc_geno_arr = chr20_data.get_lanc_geno(np.arange(10, 14, dtype=np.uint32))
    lanc_geno_true = np.asarray(
        [
            [[0, 0], [0, 0], [0, 0], [0, 0]],
            [[0, 0], [0, 0], [0, 1], [0, 0]],
            [[0, 0], [0, 0], [0, 1], [0, 0]],
            [[0, 0], [0, 0], [0, 1], [0, 0]],
            [[0, 0], [0, 1], [0, 0], [0, 0]],
            [[0, 0], [1, 0], [0, 0], [0, 0]],
            [[0, 1], [0, 0], [0, 0], [0, 0]],
            [[0, 0], [0, 0], [1, 1], [1, 0]],
            [[0, 0], [0, 0], [1, 1], [0, 0]],
            [[0, 0], [1, 0], [1, 0], [0, 0]],
            [[0, 0], [0, 0], [0, 0], [0, 0]],
            [[0, 1], [0, 0], [0, 1], [0, 0]],
            [[0, 0], [0, 0], [0, 0], [0, 0]],
            [[0, 0], [1, 0], [0, 0], [0, 0]],
            [[0, 0], [0, 0], [0, 1], [0, 0]],
            [[0, 0], [0, 0], [0, 0], [1, 0]],
            [[0, 0], [1, 0], [1, 0], [0, 0]],
            [[0, 0], [0, 0], [0, 0], [0, 0]],
            [[0, 0], [0, 1], [0, 2], [0, 0]],
            [[0, 0], [0, 0], [0, 0], [0, 0]],
        ],
        dtype=np.int32,
    )

    np.testing.assert_equal(lanc_geno_arr, lanc_geno_true)
