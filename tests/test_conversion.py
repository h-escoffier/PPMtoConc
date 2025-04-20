from ppm_to_conc import conversion


def test_conv_ppm_to_g():
    ppm = 459 
    molecular_weight = 106874  # Dalton

    result = conversion.ppm_to_grams(ppm, molecular_weight)

    # Check that the result is numeric and positive
    assert isinstance(result, float)
    assert result > 0


def test_calculate_sequence_mw():
    sequence_test = 'MAAVAAVAAR'
    mw = conversion.calculate_sequence_mw(sequence_test)

    # Check that the result is numeric and positive
    assert isinstance(mw, float)
    assert mw > 0
