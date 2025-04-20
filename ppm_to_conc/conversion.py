from Bio.SeqUtils import molecular_weight
from Bio.Seq import Seq


def calculate_sequence_mw(sequence):
    return molecular_weight(Seq(sequence), seq_type='protein')


def ppm_to_grams(ppm, molecular_weight_dalton):
    relative_abundance = ppm / 1e6
    return (relative_abundance / 6.022e23) * molecular_weight_dalton
