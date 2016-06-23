#!/usr/bin/env python

import os
import sys
import time
import logging
import argparse
from subprocess import call
from collections import defaultdict

fnull = open(os.devnull,'w')

prehelp='''
NAME
    filter_sequence_by_RDP

VERSION
    -

SYNOPSIS
    filter_sequence_by_RDP [options] -f gene_assembly 

DESCRIPTION
    It is to filter sequences with low scores for the
    assignment at a specified taxonomic rank.
'''

posthelp='''
OUTPUT
    A filtered FASTA file.

EXAMPLES
    -
 
SEE ALSO
    -

AUTHORS
    Feng Zeng
'''


phylogeny_level = ['domain','phylum','class','order','family','genus']
phylogeny_level_prefix = {'domain':'0__',\
                          'phylum':'1__',\
                          'class':'2__',\
                          'order':'3__',\
                          'family':'4__',\
                          'genus':'5__'}
phylogeny_rank_index = {'domain':0,'d':0,\
                        'phylum':1,'p':1,\
                        'class':2,'c':2,\
                        'order':3,'o':3,\
                        'family':4,'f':4,\
                        'genus':5,'g':5}

def parse_rdp_result(rdp_result):
    levels = phylogeny_level
    gene_lineage = defaultdict(list)
    with open(rdp_result) as f:
        for line in f:
            fields = line.strip().split('\t')
            gene = fields[0]
            taxon,level,score = '','',0
            lineage = []
            terms = []
            for t in fields[1:]:
                if t in levels:
                    level = t
                elif len(level)==0:
                    terms += [t]
                else:
                    score = float(t)
                    taxon = ' '.join([x.replace('"','') for x in terms])
                    lineage += [(phylogeny_level_prefix[level]+taxon,score)]
                    taxon,level,score = '','',0
                    terms = []
            gene_lineage.setdefault(gene,lineage)
    return gene_lineage
    
def run_rdp_classifier(rdp_classifier,gene_seq):
    prefix = os.path.splitext(gene_seq)[0]
    rdp_out = prefix+'_rdp_result.txt'
    cmd = ['java','-jar',rdp_classifier,'-f','fixrank','-o',rdp_out,gene_seq]
    call(cmd,stderr=fnull)
    return rdp_out

def main():
    global opts
    # run rdp classifier
    rdp_results = run_rdp_classifier(opts.rdp_classifier,opts.fasta)

    # parse rdp results
    gene_lineage = parse_rdp_result(rdp_results)

    # rank index
    rank_idx = phylogeny_rank_index[opts.rank]

    # filtration
    out = os.path.splitext(opts.fasta)[0]+'_filtered.fasta'
    fout = open(out,'w')

    cmd = ['samtools','faidx',opts.fasta]
    for gene,lineage in gene_lineage.iteritems():
        taxon,score = lineage[rank_idx]
        if score >= opts.thresh:
            cmd += [gene]

    call(cmd,stdout=fout,stderr=fnull)
    fout.close()
    call(['samtools','faidx',out])

    call(['rm',rdp_results])
        

if __name__=='__main__':
    try:
        t_start = time.time()
        # set command-line parser
        parser = argparse.ArgumentParser(description=prehelp,epilog=posthelp,\
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument('-f','--fasta',help='gene sequences',metavar='GENE_SEQ',required=True)
        parser.add_argument('-C','--rdp-classifier',help='path to the RDP classifier',dest='rdp_classifier',\
                            default='/home/fzeng/Tool/MetagenomeTools/RDPTools/classifier.jar',type=str)
        parser.add_argument('-r','--rank',help='specific taxonomic rank where the score will be used [g]',\
                            default='g',dest='rank',type=str)
        parser.add_argument('-t','--thresh',help='specify phylogeny assignment threshold for copy number correction [0.6]',\
                            dest='thresh',default=0.6,type=float)

        # parse options
        opts = parser.parse_args()
        opts.fasta = os.path.abspath(opts.fasta)

        # check rdp classifier
        if not os.path.exists(opts.rdp_classifier):
            print >>sys.stderr,"Need to specify the RDP classifier for copy number correction!"
            print >>sys.stderr,"Exit the program!"
            sys.exit(0)

        # run the program
        main()

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
