from ppm_to_conc import conversion


def test_conv_ppm_to_g():
    # Example values
    ppm = 50000
    molecular_weight = 50000  # daltons

    result = conversion.ppm_to_grams(ppm, molecular_weight)

    # Check that the result is numeric and positive
    assert isinstance(result, float)
    assert result > 0
