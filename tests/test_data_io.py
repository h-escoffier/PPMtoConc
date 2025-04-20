import pytest
import pandas as pd
from unittest import mock
from io import StringIO

import ppm_to_conc.data_io as data_io 


def test_load_paxdb():
    mock_file = StringIO(
        "# some comment\n"
        "9606.ENSP000001\t100\n"
        "9606.ENSP000002\t200\n"
    )
    with mock.patch("builtins.open", return_value=mock_file):
        df = data_io.load_paxdb("fake_path.csv")
    assert list(df.columns) == ["string_external_id", "abundance", "ENSPID"]
    assert df.shape[0] == 2
    assert df["ENSPID"].tolist() == ["ENSP000001", "ENSP000002"]


def test_load_tsv():
    mock_file = StringIO(
        "ENSPID\tabundance\n"
        "ENSP000001\t100\n"
        "ENSP000002\t200\n"
    )
    with mock.patch("pandas.read_csv", return_value=pd.read_csv(mock_file, sep="\t")) as mock_read:
        df = data_io.load_tsv("fake.tsv")
    assert list(df.columns) == ["ENSPID", "abundance"]
    assert df.shape[0] == 2
    mock_read.assert_called_once()


def test_check_file_format_paxdb():
    mock_file = StringIO("# some comment\n9606.ENSP000001,100\n")
    with mock.patch("builtins.open", return_value=mock_file):
        assert data_io.check_file_format("fake_path.csv") == "paxdb"


def test_check_file_format_tsv():
    mock_file = StringIO("ENSPID\tabundance\nENSP000001\t100\n")
    with mock.patch("builtins.open", return_value=mock_file):
        assert data_io.check_file_format("fake_path.tsv") == "tsv"


def test_load_file_dispatch_paxdb(monkeypatch):
    monkeypatch.setattr(data_io, "check_file_format", lambda x: "paxdb")
    monkeypatch.setattr(data_io, "load_paxdb", lambda x, species=9606: "paxdb_loaded")
    result = data_io.load_file("somefile")
    assert result == "paxdb_loaded"


def test_load_file_dispatch_tsv(monkeypatch):
    monkeypatch.setattr(data_io, "check_file_format", lambda x: "tsv")
    monkeypatch.setattr(data_io, "load_tsv", lambda x: "tsv_loaded")
    result = data_io.load_file("somefile")
    assert result == "tsv_loaded"


def test_load_file_unsupported(monkeypatch):
    monkeypatch.setattr(data_io, "check_file_format", lambda x: "unknown")
    with pytest.raises(ValueError):
        data_io.load_file("somefile")


def test_get_protein_molecular_weight_fill_mean(monkeypatch):
    test_df = pd.DataFrame({"ENSPID": ["ENSP000001", "ENSP000002"]})

    # Fake mapping ENSPID -> UniProt IDs
    monkeypatch.setattr(data_io, "map_enspid_to_uniprot", lambda x: {"ENSP000001": "P1", "ENSP000002": "P2"})

    # First UniProt returns a valid weight, second returns None -> tests fill_missing
    monkeypatch.setattr(data_io, "fetch_uniprot_molecular_weight", lambda x: 42.0 if x == "P1" else None)

    # For those missing MW, simulate protein sequence fallback
    monkeypatch.setattr(data_io, "fetch_protein_sequence", lambda x: "SEQ")
    monkeypatch.setattr(data_io, "calculate_sequence_mw", lambda seq: None if seq == "SEQ" else 0)

    result = data_io.get_protein_molecular_weight(test_df, fill_missing="mean")

    assert "MolecularWeight" in result.columns
    assert result["MolecularWeight"].isnull().sum() == 0


def test_convert_ppm_dataframe():
    df = pd.DataFrame({
        "ENSPID": ["E1", "E2"],
        "UniProtID": ["P1", "P2"],
        "MolecularWeight": [50000, 60000],
        "abundance": [100, 200]
    })

    with mock.patch("ppm_to_conc.data_io.ppm_to_grams", side_effect=lambda ppm, mw: ppm * mw * 1e-9):
        converted = data_io.convert_ppm_dataframe(df)
    
    assert "Mass_g_per_gDCW" in converted.columns
    assert abs(converted["Mass_g_per_gDCW"].sum() - 0.505717472) < 1e-6  # should sum to total_prot_content
