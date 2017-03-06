# ccbb libraries
from dual_crispr_pipeliner import generate_params, run_pipeline

__author__ = 'Amanda Birmingham'
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


def run_mali_pipeline(fastq_set_name, human_readable_name, library_name, num_processors, dirs_dict):
    # *********************************************************
    # Change these lines only when adding or changing a construct library from which you're analyzing data.
    # Note that the file named here is expected to be in tab-delimited text format with no quotes around text values,
    # and that the first line will be skipped (as a header row).  For this file, provide in the order shown below the
    # 0-based column indices of each of the following required columns:
    # Sequence ID
    # Gene A gRNA Sequence
    # Gene B gRNA Sequence

    if library_name == "CV4":
        # for CV4 library:
        spacers_file_name = "CV4_2spacers.txt"
        col_indices = "1,6,10"
        min_trimmed_grna_len = 19
        max_trimmed_grna_len = 21
    elif library_name == "MV4":
        # for MV4 library:
        spacers_file_name = "Metabolism_dual_spacers.txt"
        col_indices = "1,6,10"
        min_trimmed_grna_len = 19
        max_trimmed_grna_len = 21
    elif library_name == "CRV4":
        # for CRV4 library:
        spacers_file_name = "Cancer_rep2_CRV4.txt"
        col_indices = "1,3,5"
        min_trimmed_grna_len = 19
        max_trimmed_grna_len = 24
    else:
        raise ValueError("Unrecognized library name: {0}".format(library_name))

    # *********************************************************
    # Change these lines only when changing the scaffold sequence with which your library was constructed.
    # Note that these should be the FULL scaffold sequence, even if the sequencing will not capture the whole 3'
    # scaffold sequence; the cutadapt settings used will still identify a significantly truncated 3' match.
    full_5p_r1 = "TATATATCTTGTGGAAAGGACGAAACACCG"
    full_5p_r2 = "CCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC"
    full_3p_r1 = "GTTTCAGAGCTATGCTGGAAACTGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAG"
    full_3p_r2 = "CAAACAAGGCTTTTCTCCAAGGGATATTTATAGTCTCAAAACACACAATTACTTTACAGTTAGGGTGAGTTTCCTTTTGTGCTGTTTTTTAAAATA"

    # *********************************************************
    # Decide on these values before beginning screen analysis, set them once, and then leave them alone
    # unless you change your analysis approach
    len_of_seq_to_match = 19
    num_allowed_mismatches = 1

    # *********************************************************
    # DON'T change anything below here unless you *really* know what you're doing
    shared_params = generate_params(fastq_set_name, human_readable_name, num_processors, dirs_dict, spacers_file_name,
                                    col_indices, min_trimmed_grna_len, max_trimmed_grna_len, full_5p_r1, full_5p_r2,
                                    full_3p_r1, full_3p_r2, len_of_seq_to_match, num_allowed_mismatches)

    run_pipeline(dirs_dict, shared_params)
