#!/usr/bin/env python3

STANDARD_GENETIC_CODE = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def get_reverse(sequence: str) -> str:
    """
    Return the reverse of a nucleotide sequence string (uppercased).

    Examples
    --------
    >>> get_reverse('AUGC')
    'CGUA'
    """
    return sequence[::-1].upper()


def get_complement(sequence: str) -> str:
    """
    Return the complement of a nucleotide sequence string (uppercased).

    Examples
    --------
    >>> get_complement('AUGC')
    'UACG'
    """
    complement_map = str.maketrans('AUGCaugc', 'UACGuacg')
    return sequence.translate(complement_map).upper()


def reverse_and_complement(sequence: str) -> str:
    """
    Return the reverse complement of a nucleotide sequence string.

    Examples
    --------
    >>> reverse_and_complement('AUGC')
    'GCAU'
    """
    return get_reverse(get_complement(sequence))


def translate_sequence(rna_sequence: str, genetic_code: dict = None) -> str:
    """
    Translate a sequence of RNA into amino acids, reading codons from
    position 0. Stop at the first stop codon ('*').

    Examples
    --------
    >>> translate_sequence('GUCGAA')
    'VE'
    """
    if genetic_code is None:
        genetic_code = STANDARD_GENETIC_CODE

    rna_sequence = rna_sequence.upper()
    amino_acids = []
    for i in range(0, len(rna_sequence) - 2, 3):
        codon = rna_sequence[i:i+3]
        aa = genetic_code.get(codon, '')
        if aa == '*':
            break
        amino_acids.append(aa)
    return ''.join(amino_acids)


def get_all_translations(rna_sequence: str, genetic_code: dict = None) -> list:
    """
    Return a list of amino acid sequences from all AUG-initiated reading
    frames across both the forward strand and its reverse complement.

    Returns
    -------
    list of str
        Non-empty amino acid sequences starting from each AUG found.
    """
    if genetic_code is None:
        genetic_code = STANDARD_GENETIC_CODE

    def translations_from(seq):
        found = []
        start = 0
        while True:
            pos = seq.find('AUG', start)
            if pos == -1:
                break
            peptide = translate_sequence(seq[pos:], genetic_code)
            if peptide:
                found.append(peptide)
            start = pos + 1
        return found

    seq_upper = rna_sequence.upper()
    results = translations_from(seq_upper)
    results.extend(translations_from(reverse_and_complement(seq_upper)))
    return results


def get_longest_peptide(rna_sequence: str, genetic_code: dict = None) -> str:
    """
    Return the longest peptide found across all AUG-initiated reading frames.
    Returns empty string if no peptides found.
    """
    if genetic_code is None:
        genetic_code = STANDARD_GENETIC_CODE

    peptides = get_all_translations(rna_sequence, genetic_code)
    if not peptides:
        return ''
    return max(peptides, key=len)


if __name__ == '__main__':
    genetic_code = STANDARD_GENETIC_CODE
    rna_seq = 'AUGGUACCCGGGUUUAAA'
    print(get_all_translations(rna_seq, genetic_code))
    print(get_longest_peptide(rna_seq, genetic_code))
