# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 14:11:32 2022

@author: Shenied E Maldonado Guerra
"""

import sys
import csv
import numpy as np


def read_file():
    """
    Reads input file.

    Returns
    -------
    sequences : 2D list
        List containing pairs of sequences found in rows
        on input file.

    """

    fname = sys.argv[1]

    with open(fname, 'r') as file:

        reader = csv.reader(file)

        # Skip header
        next(reader)

        # Extract sequences
        sequences = []

        try:
            for row in reader:
                sequences.append(row)
        except:
            pass

        return sequences


def match_score(A, B, match=1, missmatch=-1):
    if A == B:
        return match
    else:
        return missmatch


def nw_scoringMatrix(s1, s2, match=1, missmatch=-1, gap=-2):
    """

    Parameters
    ----------
    s1 : string
        first sequence.
    s2 : string
        second sequence.
    match : integer, optional
        score for match between two residues. The default is 1.
    missmatch : integer, optional
        score for missmatch between two residues. The default is -1.
    gap : gap penalty, optional
        score for including a gap. The default is -2.

    Returns
    -------
    None.

    """
    nx = len(s1) + 1  # rows
    ny = len(s2) + 1  # columns

    # Initialization: Scoring Matrix
    F = np.zeros((nx, ny))
    F[:, 0] = np.linspace(start=0, stop=len(s1) * -2, num=nx)
    F[0, :] = np.linspace(start=0, stop=len(s2) * -2, num=ny)

    # Computing scoring matrix
    for i in range(len(s1)):  # rows
        for j in range(len(s2)):  # columns

            # Calculate maximum score: Diagonal, Above or Left
            diag_score = F[i, j] + match_score(s1[i], s2[j])
            up_score = F[i+1, j] + gap
            left_score = F[i, j+1] + gap

            largest = max([diag_score, up_score, left_score])

            # Fill out scoring matrix
            F[i+1, j+1] = largest

    return nw_traceback(F, s1, s2)


def nw_traceback(F, s1, s2, match=1, missmatch=-1, gap=-2):
    """

    Parameters
    ----------
    F : Scoring Matrix
        Matrix or 2D list with the same amount of rows as there are letters in
        the first sequence & the same amount of columns are there are letters
        in the second sequence.
    s1 : String
        First sequence.
    s2 : String
        Second sequence.

    Returns
    -------
    results : list
        Original two sequences, their alignment & their score.

    """
    # Iterators
    i = len(s1)
    j = len(s2)

    rseq1 = ''
    rseq2 = ''

    # Final score
    score = F[i, j]

    while i > 0 and j > 0:

       # Diagonal, Up, Left
       options = [F[i-1, j-1], F[i-1, j], F[i, j-1]]

       largest = max(options)

       if largest == F[i-1, j-1]:
           # Move diagonally
           rseq1 += s1[i-1]
           rseq2 += s2[j-1]
           i -= 1
           j -= 1

       elif largest == F[i-1, j]:
           # Move Up
           rseq1 += s1[i-1]
           rseq2 += '-'
           i -= 1

       elif largest == F[i, j-1]:
           # Move Left
           rseq1 += '-'
           rseq2 += s2[j-1]
           j -= 1

    # Remaining nucleotides in either sequence
    if i > 0:
        while i > 0:
            rseq1 += s1[i-1]
            i -= 1

    if j > 0:
        while j > 0:
            rseq2 += s2[j-1]
            j -= 1

    # Reverse strings
    rseq1 = rseq1[::-1]
    rseq2 = rseq2[::-1]

    return [s1, s2, rseq1, rseq2, score]


def write_file():
    """
    
    Runs Needleman-Wunsch algorithm on input file's data and produces a 
    results.csv file.'

    Returns
    -------
    None.

    """

    f = open('results.csv', 'a', newline='')

    sequences = read_file()

    writer = csv.writer(f)

    writer.writerow(['sequence1', 'sequence2',
                    'alignment text', 'alignment score'])

    for pair in sequences:
        result = nw_scoringMatrix(pair[0], pair[1])

        writer.writerow([result[0], result[1], result[2] +
                        "\n" + result[3], result[4]])

    f.close()


write_file()
