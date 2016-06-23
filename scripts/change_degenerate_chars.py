#!/usr/bin/env python
'''
Created by Feng Zeng    Jan 23, 2016
'''

import os
import sys
import argparse
from random import choice

prehelp='''
NAME
    change_degenerate_chars

VERSION
    -

SYNOPSIS
    change_degenerate_chars FASTA_FILE

DESCRIPTION
    It is used to change degenerate chars to non-degenerate chars 
    for the input sequence(s).
'''

posthelp='''
OUTPUT
    The non-degenerate character sequence(s) will be output to
    the stdout stream.

EXAMPLES
    -

SEE ALSO
    -

AUTHORS
    Feng Zeng
'''


IUPAC_CODE = {"-": "-",
              ".": "-",
              "A": "A",
              "C": "C",
              "G": "G",
              "T": "T",
              "U": "T",
              "R": "AG", #        Purine (A or G)
              "Y": "CT", #        Pyrimidine (C, T, or U)
              "M": "CA", #        C or A
              "K": "TG", #        T, U, or G
              "W": "TA", #        T, U, or A
              "S": "CG", #        C or G
              "B": "CTG",#        C, T, U, or G (not A)
              "D": "ATG",#        A, T, U, or G (not C)
              "H": "ATC",#        A, T, U, or C (not G)
              "V": "ACG",#        A, C, or G (not T, not U)
              "N": "ACTG" }#      Any base (A, C, G, T, or U)

iupac_map = IUPAC_CODE.get


def change_degenerate_chars(fasta_file):
    '''
    It is the function for mapping degenerate chars to
    non-degenerate chars.
    '''

    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                print line.strip()
            else:
                new_line = ''.join([choice(iupac_map(x.upper(),'ACGT')) for x in line.strip()])
                print new_line


if __name__=='__main__':
    parser = argparse.ArgumentParser(usage='',description=prehelp,epilog=posthelp,\
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('fasta',help='Input sequences, FASTA format',metavar='FASTA_FILE')
    
    # parse options
    opts = parser.parse_args()

    # execute
    change_degenerate_chars(opts.fasta)
