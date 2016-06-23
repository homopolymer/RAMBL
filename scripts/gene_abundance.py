#!/home/fzeng/anaconda/bin/python
"""
Calculate gene abundance.

Change Log
==========
Sep 24, 2015    Feng Zeng    Add coverage ratio and print out
"""

from __future__ import division
import os
import sys
import csv
import logging
import argparse
import numpy as np
import time
from collections import defaultdict
import textwrap,traceback
from itertools import izip

def cal_gene_abundance(breadth,depth):
    """abundance equation: sum breadthxdepth over all segments
    
    Examples
    --------
    1. A gene is divided into 3 segments, the corresponding breadths
       of which are 0.1,0.6,0.3.  The depth of these three segments 
       are 0,8,10.  Therefore, the abundance is calculated as

           abundance = (0.1)^1.5x0 + (0.6)^1.5x8 + (0.3)^1.5x10 

    """
    lam = 1.0
    return np.dot(np.power(breadth,lam),depth)

def cal_gene_coverage(breadth,depth):
    """calculate the coverage ratio of a gene"""
    ratio = []
    for r,a in izip(breadth,depth):
        if a>0:
            ratio += [r]
    return sum(ratio)


def main():
    """the main interface"""
    global opts

    # TODO: load gene info
    if opts.verbose:
        logging.info('load gene info from %s' % os.path.abspath(opts.gi[0]))
    gene = defaultdict(int)
    with open(opts.gi[0],'r') as f:
        r = csv.reader(f,delimiter="\t")
        for row in r:
            gene[row[0]] = int(row[1])

    # TODO: calculate gene abundance
    if opts.verbose:
        logging.info('process records in %s' % os.path.abspath(opts.cov[0]))
    with open(opts.cov[0],'r') as f:
        gn = None             # gene name
        breadth = list()      # segment breadth
        depth = list()        # segment depth
        r = csv.reader(f,delimiter="\t")
        for row in r:
            if gn is not None and gn != row[0]:
                # calculate gene abundance and report
                abun = cal_gene_abundance(breadth,depth)
                ratio = cal_gene_coverage(breadth,depth)
                print "%s\t%d\t%d\t%f\t%f" % (gn,1,gene[gn],abun,ratio)
                # clean state
                breadth = list()
                depth = list()
            # record state
            gn = row[0]
            breadth.append((float(row[2])-float(row[1])+1.)/gene[gn])
            depth.append(float(row[3]))
        # last record
        abun = cal_gene_abundance(breadth,depth)
        ratio = cal_gene_coverage(breadth,depth)
        print "%s\t%d\t%d\t%f\t%f" % (gn,1,gene[gn],abun,ratio)
    
    if opts.verbose:
        logging.info('done')
            

if __name__ == "__main__":
    try:
        t_start = time.time()
        # TODO: set command-line parser
        parser = argparse.ArgumentParser(description=globals()['__doc__'],epilog=textwrap.dedent('''\
                                         Examples\n\
                                         --------\n\
                                         1. gene_abundance to_gg_99_otus.cov gg_99_otus.fasta.fai\n \
                                         '''),\
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument('cov',help='gene depth file',nargs=1,metavar='DEPTH')
        parser.add_argument('gi',help="gene index file",nargs=1,metavar="GENE")
        parser.add_argument('-v',dest='verbose',action='store_true',default=False,help='verbose output')
        opts = parser.parse_args()

        # TODO: set logging format
        logging.basicConfig(format='[%(asctime)s] %(levelname)s : %(message)s', level=logging.INFO)   

        # TODO: goto main interface
        main()

        # TODO: complete
        if opts.verbose:
            logging.info('elapsed time is {} minutes'.format(round((time.time()-t_start)/60.,5)))

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
