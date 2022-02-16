"""
Library for finding potential genes in a strand of DNA.
"""
from tqdm import tqdm_notebook as tqdm
from helpers import amino_acid, shuffle, load_fasta_file
tqdm().pandas()

def get_complement(nucleotide):
    """
    Takes a string nucleotide consisting of a single character
    A, T, C, or G, each representing the DNA nucleotides
    adenine, thymine, cytosine, and guanine,
    and return a string (also consisting of a single character)
    representing the complementary nucleotide.
    In DNA, A and T are complementary to each other,
    as are C and G

    Args:
    - nucleotide: a string that represents the DNA nucleotides
    in a certain sequence with A, T, C, or G.

    Returns:
    - a string (consisting of a single character)
    representing the complementary nucleotide.
    """
    complement = ""
    if nucleotide in 'A':
        complement = 'T'
    elif nucleotide in 'T':
        complement = 'A'
    elif nucleotide in 'G':
        complement = 'C'
    else:
        complement = 'G'
    return complement

def get_reverse_complement(strand):
    """
    Takes a string strand representing a single strand of DNA
    (consisting entirely of the characters A, T, C, and G)
    and returns a string representing a complementary strand of DNA

    Args:
    - strand: a string that represents a string with multiple
    characters that are a single strand of DNA

    Returns:
    - a string representing a complementary strand of DNA
    """
    reverse_sequence = ''

    for base in reversed(strand):
        reverse_sequence += get_complement(base)
    return reverse_sequence

def rest_of_orf(strand):
    """
    Takes a string strand representing a strand of DNA
    that begins with a start codon and returns the sequence
    of nucleotides representing the rest of the ORF,
    up to but not including the next stop codon
    which include 'TAA', 'TAG', 'TGA.'

    If there is no stop codon, the function
    returns the entire strand.

    Args:
    - strand: a string representing a strand of DNA

    Returns:
    -  a string that begins with a start codon
    and returns the rest of the opening reading frame
    but not including the next stop codon.
    """
    # Makes the strand list which is the strand
    # into a list for each triplet
    strand_list = [strand[nucleotide:nucleotide+3]
                   for nucleotide in range(0, len(strand), 3)]

    # For each value in the strand list
    for codon_index, codon in enumerate(strand_list):
        # if the return value of function "amino_acid"
        # is found in stop_codons list, then stop loop
        if amino_acid(codon) == "*":
            strand_list = strand_list[0:codon_index]
            break

    # Return the list as one big string
    return ''.join(strand_list)


def find_all_orfs_one_frame(strand):
    """
    Takes a string strand representing a strand of DNA
    and returns a list of strings representing
    all in-frame ORFs found in that strand

    Args:
    - strand: a string representing a strand of DNA

    Returns:
    - a list of strings representing
    all in-frame ORFs found in that strand.
    """
    # Make a list that will contain all orfs
    frame = ""
    in_frame_orfs = []

    # while the strand still can return an orf
    while len(strand) >= 3:
        current_codon = strand[:3]

        if current_codon != "ATG":
            strand = strand[3:]
            continue

        # run frame of orf onto strand and get
        frame = rest_of_orf(strand)

        strand = strand[len(frame):]

        # append the strand to the string when completed
        in_frame_orfs.append(frame)

    # return string when completed
    return in_frame_orfs

def find_all_orfs(strand):
    """
    Takes a string strand representing a strand of DNA
    and returns a list of strings representing
    not only in-frame ORFs, but also ORFs
    found one or two nucleotides from the start of strand

    Args: a string "strand" representing a strand of DNA

    Returns: a list of strings representing
    all ORFs found in that strand
    """
    all_orfs = []

    # go through every codon from the start
    # and find all the orfs in that frame
    # when completed, shift the strand by one
    for start in range(3):
        all_orfs += find_all_orfs_one_frame(strand[start:])

    return all_orfs

def find_all_orfs_both_strands(strand):
    """
    Takes a string strand representing a strand of DNA
    and returns a list of strings representing all ORFs
    found in strand or its reverse complement.

    Args: a string strand representing a strand of DNA

    Returns: a list of strings representing all ORFs
    found in strand or its reverse complement.
    """
    # create an empty list that will contain all the ORFs
    all_frame_orfs_both_strands = []

    # run the function twice for the original and reverse strand
    for _ in range(2):
        # add the function output into the list
        all_frame_orfs_both_strands.append(find_all_orfs(strand))
        strand = get_reverse_complement(strand)

    # combine the nested list into a single list
    all_frame_orfs_both_strands = [
        orf for orf_list in all_frame_orfs_both_strands for orf in orf_list]

    return all_frame_orfs_both_strands


def find_longest_orf(strand):
    """
    Finds the longest ORF found in either the strand
    or its reverse complement.

    Args: a string strand representing a strand of DNA

    Returns: the length of the longest ORF found in either
    that strand or its reverse complement.
    """
    # run find_all_orfs_both_strands into the list_of_strings
    list_of_orfs = find_all_orfs_both_strands(strand)
    longest_orf_length = 0
    longest_orf = ""
    # Go through the whole list
    for orf in list_of_orfs:
        # Checking for the longest word(string)
        if len(orf) > longest_orf_length:
            longest_orf_length = len(orf)
            longest_orf = orf
    return longest_orf

def noncoding_orf_threshold(strand, num_trials):
    """
    Takes a string strand representing a strand of DNA
    and a positive integer num_trials representing a number
    of trials to run. For each of these trials, it randomly
    shuffles the nucleotides in strand and finds the longest ORF
    in either the shuffled strand or its reverse complement.
    It then keeps track of the minimum length of this value over
    all of the trials and returns an integer representing
    this minimum length.

    Args:
    - a string strand representing a strand of DNA
    - a positive integer num trials representing a number
    of trials to run

    Returns: an integer representing the minimum length
    of the value over all of the trials
    """
    orfs_lengths = []

    for _ in tqdm(range(num_trials)):
        # shuffles the strand
        shuffled_strand = shuffle(strand)

        # finds the longest ORF in either the shuffled strand or its reverse complement.
        longest_orf = find_longest_orf(shuffled_strand)

        # keeps track of the minimum length of this value over all of the trials
        length_of_longest_orf = len(longest_orf)

        # returns an integer representing this minimum length
        orfs_lengths.append(length_of_longest_orf)

    # returns an integer representing this minimum length
    return min(orfs_lengths)


def encode_amino_acids(orf):
    """
    Determines the sequence of amino acids these ORFs encode
    Each of these sequences is potentially a protein of interest
    in the genome, and will be a candidate for further analysis.

    It is possible that an ORF's length will not be an exact multiple of 3,
    due to reading an ORF from an offset.
    In this case, there will be 1 or 2 nucleotides left at the end
    your implementation should simply ignore these.

    Args: a string orf representing a strand of DNA that is an ORF

    Returns: a string representing the sequence of amino acids
    with each amino acid written as its one-letter symbol.
    """
    # empty string for the amino acids
    amino_acids = ""
    # run the amino acid function for each codon
    # and add it to the string
    for triple in range(0, len(orf), 3):
        codon = orf[triple:triple + 3]
        amino_acids += amino_acid(codon)

    # return the string when completed
    return amino_acids


def find_genes(path):
    """
    Finds potential protein-coding genes by

    1. Loading the string of nucleotides from the file.
    2. Determining the threshold length to use as a cutoff
    for coding ORFs, using 1,500 trials of shuffling
    the nucleotide sequence.
    3. Find all ORFs in both strands of the sequence.
    4. For any ORF longer than the cutoff length,
    translate the ORF to a sequence of amino acids.

    Args: the location of a file in a FASTA format

    Returns: list of all such amino acid sequences.
    """
    # to be returned
    amino_acid_sequences = []

    # Loading the string of nucleotides from the file.
    genes = load_fasta_file(path)

    # Determining the threshold length to use as a cutoff
    # for coding ORFs, using 1,500 trials of shuffling
    # the nucleotide sequence.
    threshold_length = noncoding_orf_threshold(genes, 1500)

    # Find all ORFs in both strands of the sequence.
    all_orfs_list = find_all_orfs_both_strands(genes)

    # For any ORF longer than the cutoff length,
    # translate the ORF to a sequence of amino acids.
    # tqdm is used for the progress bar in the
    # 'sars-cov-2' file
    for orf in tqdm(all_orfs_list):
        if len(orf) > threshold_length:
            protein_orf = encode_amino_acids(orf)
            amino_acid_sequences.append(protein_orf)

    return amino_acid_sequences
