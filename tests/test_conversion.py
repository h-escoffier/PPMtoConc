import ppm_to_conc.conversion as conversion


def test_calculate_sequence_mw():
    sequence_test = 'MAAVAAVAAR'
    mw = conversion.calculate_sequence_mw(sequence_test)

    assert isinstance(mw, float)
    assert mw > 0


def test_ppm_to_grams():
    ppm = 500 
    molecular_weight = 100000 

    result = conversion.ppm_to_grams(ppm, molecular_weight)

    assert isinstance(result, float)
    assert result > 0
