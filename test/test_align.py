# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    nw = NeedlemanWunsch("substitution_matrices/BLOSUM62.mat", -10, -1)
    score, align1, align2 = nw.align(seq1, seq2)

    #set the manual alignemnt to compare to
    manual_align1 = "MYQR"
    manual_align2 = "M-QR"
    print (manual_align1, manual_align2)
    print(align1, align2)
    #print the matrices
    print(nw._align_matrix)
    print(nw._gapA_matrix)
    print(nw._gapB_matrix)
    print(nw._back)
    #assert that the alignment is correct
    assert align1 == manual_align1
    assert align2 == manual_align2
    print(score)

    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")

    nw = NeedlemanWunsch("substitution_matrices/BLOSUM62.mat", -10, -1)
    score, align1, align2 = nw.align(seq3, seq4)
    manual_align1 = "MAVHQLIRRP"
    manual_align2 = "M---QLIRHP"
    print (manual_align1, manual_align2)
    print(align1, align2)
    #assert that the alignment and backtracing are correct
    assert align1 == manual_align1
    assert align2 == manual_align2
    #assert that the score is correct
    assert score == 17
    print(score)




