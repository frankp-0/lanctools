import pytest
import numpy as np
import pandas as pd
from lanctools import LancData
from lanctools import convert_to_lanc
from lanctools.core import _parse_lanc_line, _read_lanc


@pytest.fixture
def chr20_data():
    dataset = LancData(
        plink_prefix="tests/data/chr20", lanc_file="tests/data/chr20.lanc"
    )
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


def test_convert_flare(tmp_path):
    tmp_lanc_path = tmp_path / "test_flare.lanc"
    convert_to_lanc(
        file="tests/data/chr20.flare.anc.vcf.gz",
        file_fmt="FLARE",
        plink_prefix="tests/data/chr20",
        output=tmp_lanc_path,
    )

    with (
        open("tests/data/chr20.lanc", "r", encoding="utf-8") as true_lanc,
        open(tmp_lanc_path, "r", encoding="utf-8") as test_lanc,
    ):
        assert true_lanc.read() == test_lanc.read()


def test_get_info(chr20_data):
    nvar = chr20_data.pvar.get_variant_ct()
    df_info = chr20_data.get_info(np.arange(nvar))
    df_true = pd.read_json(
        "tests/data/chr20_info.json",
        dtype={"chrom": str, "pos": np.uint32, "ref": str, "alt": str, "rsid": str},
    )
    pd.testing.assert_frame_equal(df_info, df_true)


# TODO: test out of bounds indices
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
