#!/usr/bin/env python
'''

Change Log
==========
Dec 8, 2015    Feng Zeng    Create it
'''

import os
import sys
import csv
import time
import logging
import argparse
import subprocess
from collections import defaultdict
from skbio import SequenceCollection,DNA,Sequence
from random import choice

prehelp = """
NAME
    build_silva_data

VERSION
    -

SYNOPSIS
    build_silva_data [options] <SILVA_SEQ> <SILVA_ALIGN>

DESCRIPTION
    It is to canonize sequence id and compute the phylogeny
    tree for SILVA, and prepare SILVA data for RBRA. 
"""

posthelp = """
OUTPUT
    A new fasta file wherein sequence tags are canonical.
    A tree file corresponding to the SILVA alignment.

EXAMPLES
    -

SEE ALSO
    -

AUTHORS
    Feng Zeng 
"""

replace_dict = {"-": "-",
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

def proc_seq(silva_seq):
    '''
    load SILVA sequence and change accession id
    '''
    global opts

    if opts.verbose:
        logging.info('processing SILVA sequence')

    get = replace_dict.get

    outfile = '%s.fasta'%opts.prefix
    oh = open(outfile,'w')

    i = 0
    with open(silva_seq) as infile:
        for line in infile:
            line = line.rstrip()
            if line.startswith('>'):
                if i%5000==0 and opts.verbose:
                    print >>sys.stderr,i
                i += 1
                fields = line.split()
                line = fields[0].replace('.','_')
            else:
                line = ''.join([choice(get(x,'ACTG')) for x in line])
            oh.write('%s\n'%line)
  
    oh.close()
    

def proc_align(silva_align):
    '''
    load SILVA alignment and change accession id
    '''
    global opts

    if opts.verbose:
        logging.info('processing SILVA alignment')

    outfile = '%s_align.fasta'%opts.prefix
    oh = open(outfile,'w')

    get = replace_dict.get

    i = 0
    with open(silva_align) as infile:
        for line in infile:
            line = line.rstrip()
            if line.startswith('>'):
                if i%5000==0 and opts.verbose:
                    print >>sys.stderr,i
                i += 1
                fields = line.split()
                line = fields[0].replace('.','_')
            else:
                line = ''.join([choice(get(x,'ACTG')) for x in line])
            oh.write('%s\n'%line)

    oh.close()
  


def build_tree():
    '''
    build phylogeny tree by FastTree
    '''
    global opts
    
    if opts.verbose:
        logging.info('use FastTree to compute phylogeny tree')

    os.environ["OMP_NUM_THREADS"] = str(opts.cores)
    cmd = ['FastTreeMP','-nt','%s_align.fasta'%opts.prefix]
    with open('%s.tre'%opts.prefix,'w') as fh:
        subprocess.call(cmd,stdout=fh)



if __name__ == "__main__":
    try:
        t_start = time.time()
        # set command-line parser
        parser = argparse.ArgumentParser(usage='',description=prehelp,epilog=posthelp,\
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument('seq',help='SILVA sequences',metavar='SILVA_SEQ')
        parser.add_argument('align',help='SILVA alignment',metavar='SILVA_ALIGN')
        parser.add_argument('-p','--prefix',help='output filename prefix [silva_data]',\
                            default='silva_data',dest='prefix',metavar='STR')
        parser.add_argument('-s','--start-from',help='specify the starting step [0]',default=0,type=int,dest='start')
        parser.add_argument('-c','--cores',help='specify the CPU cores [1]',default=1,type=int,dest='cores')
        parser.add_argument('-v','--verbose',help='verbose output',dest='verbose',action='store_true',default=False)

        # parse options
        opts = parser.parse_args()

        # set logging
        logging.basicConfig(format='[%(asctime)s] %(levelname)s : %(message)s', level=logging.INFO)

        # 1. process sequences
        if opts.start<=0:
            proc_seq(opts.seq)

        # 2. process alignments
        if opts.start<=1:
            proc_align(opts.align)

        # 3. build tree
        if opts.start<=2:
            build_tree()
        
        # complete
        if opts.verbose:
            logging.info('elapsed time is %.5f minutes' % ((time.time()-t_start)/60.))

        sys.exit(0)
    except KeyboardInterrupt,e:
        raise e
    except SystemExit,e:
        raise e
    except Exception,e:
        logging.exception('ERROR, UNEXPECTED EXCEPTION')
        logging.exception(e)
        traceback.print_exc()
        os._exit(1)
