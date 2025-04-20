from ppm_to_conc import conversion


def test_calculate_sequence_mw():
    sequence_test = 'MAAVAAVAAR'
    mw = conversion.calculate_sequence_mw(sequence_test)

    assert isinstance(mw, float)
    assert mw > 0


def test_ppm_to_grams():
    ppm = 459 
    molecular_weight = 106874  # Dalton

    result = conversion.ppm_to_grams(ppm, molecular_weight)

    assert isinstance(result, float)
    assert result > 0
