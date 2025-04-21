import pytest
from unittest.mock import patch, MagicMock
import ppm_to_conc.api as api


@pytest.fixture
def mock_uniprot_mapping_result():
    return {
        'results': [
            {'from': 'ENSP00000354587', 'to': {'primaryAccession': 'P12345'}},
            {'from': 'ENSP00000419060', 'to': {'primaryAccession': 'Q67890'}}
        ]
    }


@patch('ppm_to_conc.api.UniProt')
def test_map_enspid_to_uniprot(mock_uniprot_class, mock_uniprot_mapping_result):
    mock_instance = mock_uniprot_class.return_value
    mock_instance.mapping.return_value = mock_uniprot_mapping_result

    test_enspids = ['ENSP00000354587', 'ENSP00000419060']
    result = api.map_enspid_to_uniprot(test_enspids, batch_size=2, sleep_time=0)

    expected = {
        'ENSP00000354587': 'P12345',
        'ENSP00000419060': 'Q67890'
    }
    assert result == expected


@patch('ppm_to_conc.api.requests.get')
def test_fetch_uniprot_molecular_weight(mock_get):
    mock_response = MagicMock()
    mock_response.ok = True
    mock_response.json.return_value = {
        'sequence': {'molWeight': 51234}
    }
    mock_get.return_value = mock_response

    uniprot_id = 'P12345'
    weight = api.fetch_uniprot_molecular_weight(uniprot_id)

    assert weight == 51234
    mock_get.assert_called_once_with(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json")


@patch('ppm_to_conc.api.requests.get')
def test_fetch_protein_sequence(mock_get):
    mock_response = MagicMock()
    mock_response.json.return_value = {'seq': 'MEEPQSDPSV'}
    mock_get.return_value = mock_response

    ensp_id = 'ENSP00000354587'
    sequence = api.fetch_protein_sequence(ensp_id)

    assert sequence == 'MEEPQSDPSV'
    expected_url = f"https://rest.ensembl.org/sequence/id/{ensp_id}?type=protein"
    mock_get.assert_called_once_with(expected_url, headers={"Content-Type": "application/json"})
