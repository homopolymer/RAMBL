#!/usr/bin/env python

import argparse
import sys,re,time
from Bio import Entrez,SeqIO
from urllib2 import HTTPError

Entrez.email ="eigtw59tyjrt403@gmail.com"

def load_fasta_id(filename):
    ids=[]
    with open(filename) as f:
        for line in f:
            if line.startswith('#') or len(line.strip())==0:
                continue
            ids.append(line.strip().split('.')[0])
    return ids


def remote_get_fasta_seqs(ids,verbose=False):
    handle = Entrez.read(Entrez.esearch(db="nucleotide", term=' '.join(ids), retmode="xml"))
    gids = handle['IdList']
    request = Entrez.epost("nucleotide",id=",".join(gids))
    result = Entrez.read(request)
    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    handle = Entrez.efetch(db="nucleotide", webenv=webEnv, query_key=queryKey, rettype="fasta", retmode="text")

    seqs = []
    i = 0
    for r in SeqIO.parse(handle,'fasta'):
        if verbose:
            print >>sys.stderr,gids[i],r.description,r.seq[:5]
        r.id = gids[i]
        seqs += [r]
        i += 1

    SeqIO.write(seqs,sys.stdout,'fasta')



def get_ncbi_fasta(filename,verbose=False):
    ids = load_fasta_id(filename)
    remote_get_fasta_seqs(ids,verbose)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('id_list_file',help='a file listing seq ids')
    parser.add_argument('-v','--verbose',dest='verbose',action='store_true',default=False,help='verbose output')
    
    opts = parser.parse_args()

    get_ncbi_fasta(opts.id_list_file,opts.verbose)
