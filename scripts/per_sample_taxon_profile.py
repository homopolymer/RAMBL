#!/usr/bin/env python
import os
import sys
import csv
import argparse
import logging
import time
import shutil
from collections import defaultdict
from subprocess import call,Popen,PIPE
from multiprocessing import cpu_count

write_devnull = open(os.devnull,'w')

prehelp='''
NAME
    sample_taxonomy_profile

VERSION
    -

SYNOPSIS
    sample_taxonomy_profile [options] gene_assembly sample_info

DESCRIPTION
    It is to invoke the RDP classifier to identify microbial taxa
    using the 16S gene sequences.
'''

posthelp='''
OUTPUT
    A tab-delimited file, of which columns are samples and rows are taxa.

EXAMPLES
    -
 
SEE ALSO
    -

AUTHORS
    Feng Zeng
'''


phylogeny_level = ['domain','phylum','class','order','family','genus']

def rdp_classify(rdp_classifier,fasta):
    cmd = ['java','-jar',rdp_classifier]
    cmd += ['-f','fixrank']
    cmd += ['-o','rdp_results']
    cmd += [fasta]
    call(cmd)

def parse_rdp_results(rdp_results,taxon_rank,thresh):
    levels = phylogeny_level
    taxa_genes = defaultdict(list)
    taxa_lineage = defaultdict()
    with open(rdp_results) as f:
        for line in f:
            fields = line.strip().split('\t')
            gene = fields[0]
            taxon,level,score = '','',0
            lineage = defaultdict()
            terms = []
            for t in fields[1:]:
                if t in levels:
                    level = t
                elif len(level)==0:
                    terms += [t]
                else:
                    score = float(t)
                    taxon = ' '.join([x.replace('"','') for x in terms])
                    taxon = taxon.strip()
                    if level==taxon_rank and score>=thresh:
                        taxa_genes.setdefault(taxon,[]).append(gene)
                    lineage.setdefault(level,taxon)
                    taxon,level,score = '','',0
                    terms = []
            taxa_lineage.setdefault(lineage[taxon_rank],lineage)
    return taxa_genes,taxa_lineage


def get_taxa_count(taxa_genes,gene_count):
    taxa_count = defaultdict(float)

    # scan through taxon
    for taxon, genes in taxa_genes.iteritems():
        count = sum([gene_count[g] for g in genes if g in gene_count])
        taxa_count.setdefault(taxon,count)

    return taxa_count
        
def get_taxa_length(taxa_genes, gene_seqs):
    if not os.path.exists(gene_seqs+'.fai'):
        call(['samtools','faidx',gene_seqs])

    gene_length = defaultdict(int)
    with open(gene_seqs+'.fai') as f:
        for line in f:
            fields = line.split()
            gene_length.setdefault(fields[0],int(fields[1]))

    taxa_length = defaultdict(int)
    for taxon,genes in taxa_genes.iteritems():
        taxa_length.setdefault(taxon,max([gene_length[g] for g in genes]))

    return taxa_length


from numpy import where,asarray

def get_taxa_copy_number(taxa, lineage, level):
    taxa_copy_number = defaultdict(float)
    copy_number_file = os.path.join(os.path.dirname(__file__),'rrnDB_RDP.tsv')
    with open(copy_number_file) as f:
        reader = csv.DictReader(f,delimiter='\t')
        for row in reader:
            taxa_copy_number.setdefault(row['name'],float(row['mean']))
    
    i  = where(asarray(phylogeny_level)==level)[0][0]

    copy_number = defaultdict(float)
    for taxon in taxa:
        if taxon in taxa_copy_number:
            copy_number.setdefault(taxon,taxa_copy_number[taxon])
        else:
            for high_level in phylogeny_level[:i][::-1]:
                high_taxon = lineage[taxon][high_level]
                if high_taxon in taxa_copy_number:
                    copy_number.setdefault(taxon,taxa_copy_number[high_taxon])
                    break

    return copy_number

def main():
    global opts
    CWD = os.getcwd()

    # make a temporary working directory
    TmpDir = os.path.join(CWD,'5_sample_taxonomy_profile')
    if os.path.exists(TmpDir):
        shutil.rmtree(TmpDir)

    os.mkdir(TmpDir)
    os.chdir(TmpDir)

    # 1 use RDP to classify
    if opts.verbose:
        logging.info('use RDP classifier to identify microbes')
    rdp_classify(opts.rdp_classifier,opts.fasta)
    taxa_genes,taxa_lineage = parse_rdp_results(os.path.join(TmpDir,'rdp_results'),opts.taxon_rank,opts.thresh)

    cmd = ['samtools','faidx',opts.fasta] 
    for genes in taxa_genes.values():
        cmd += genes
    with open(os.path.join(TmpDir,'qualified_genes.fasta'),'w') as f:
        call(cmd,stdout=f,stderr=write_devnull)

    taxa_length = get_taxa_length(taxa_genes,os.path.join(TmpDir,'qualified_genes.fasta'))

    # 3 compute gene read counts
    if opts.verbose:
        logging.info('compute read counts for genes')
    cmd = ['python',os.path.join(os.path.dirname(__file__),'per_sample_gene_profile.py')]
    cmd += ['-c',str(opts.cores)]
    cmd += ['-C',opts.rdp_classifier,'-t',str(opts.thresh)]
    cmd += ['-n']
    if opts.verbose:
        cmd += ['-v']
    cmd += [os.path.join(TmpDir,'qualified_genes.fasta'),opts.bam,opts.sample]
    call(cmd)
    
    gene_count = {}
    with open(os.path.join(TmpDir,opts.sample+'_gene_count.tsv')) as f:
        for line in f:
            if line.startswith('sample'):
                continue
            gene,count = line.strip().split()
            gene_count[gene] = float(count)
    
    # 4 taxa read count
    if opts.verbose:
        logging.info('compute sample-wise taxa read count')
    taxa_count = get_taxa_count(taxa_genes,gene_count)

    # 5 copy number correction
    taxa_copy_number = get_taxa_copy_number(taxa_genes.keys(),taxa_lineage,opts.taxon_rank)

    for taxon in taxa_count.keys():
        taxa_count[taxon] = taxa_count[taxon]/(taxa_copy_number[taxon]*taxa_length[taxon])

    # 6 normalization
    Z = sum(taxa_count.values())+1e-10
    for taxon in taxa_count.keys():
        taxa_count[taxon] = taxa_count[taxon]/Z

    fieldnames = ['sample'] + [opts.sample]
    # TODO write to disk file
    with open(os.path.join(CWD,opts.sample+'_taxa_count.tsv'),'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames,delimiter='\t')
        writer.writeheader()
        for t,c in sorted([(t,c) for t,c in taxa_count.items()]):
            writer.writerow({'sample':t,opts.sample:c})

    # back to parent directory
    os.chdir(CWD)
    # remove temporary directory
    shutil.rmtree(TmpDir)

    return None

if __name__ == "__main__":
    try:
        t_start = time.time()
        # set command-line parser
        parser = argparse.ArgumentParser(description=prehelp,epilog=posthelp,\
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument('fasta',help='gene sequences',metavar='GENE_SEQ')
        parser.add_argument('bam',help='sample reads',metavar='SAMPLE_READ')
        parser.add_argument('sample',help='sample name',metavar='SAMPLE_NAME')
        parser.add_argument('-r','--rdp-classifier',help='path to RDP classifier',\
                            default='/home/fzeng/Tool/MetagenomeTools/RDPTools/classifier.jar',\
                            dest='rdp_classifier',metavar='STR')
        parser.add_argument('-t','--thresh',help='threshold for RDP boostrap score [0.6]',\
                            default=0.6,type=float,dest='thresh',metavar='FLT')
        parser.add_argument('-L','--tax-rank',help='specify the taxonomic rank for profiling, phylum|genus [genus]',\
                            default='genus',dest='taxon_rank',metavar='STR')
        parser.add_argument('-c','--cores',help='specify the CPU cores [1]',default=1,type=int,dest='cores')
        parser.add_argument('-v','--verbose',help='verbose output',dest='verbose',action='store_true',default=False)

        # parse options
        opts = parser.parse_args()
        opts.fasta = os.path.abspath(opts.fasta)
        opts.bam = os.path.abspath(opts.bam)

        opts.cores = min([opts.cores,cpu_count()])

        # set logging
        logging.basicConfig(format='[%(asctime)s] %(levelname)s : %(message)s', level=logging.INFO)

        # check RDP classifier
        if os.path.exists(opts.rdp_classifier):
            main()

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
