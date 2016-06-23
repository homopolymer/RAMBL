#!/usr/bin/env python

import os
import sys
import time
import glob
import logging
import argparse
from subprocess import call
import extract_reads
from Bio.Blast import NCBIXML
from sets import Set
from collections import defaultdict
import csv
from multiprocessing import cpu_count,Pool
import random
from random import choice
import numpy as np
import shutil

fnull = open(os.devnull,'w')

rrnDB=os.path.join(os.path.dirname(__file__),'rrnDB_RDP.tsv')

prehelp='''
NAME
    per_sample_gene_profile

VERSION
    -

SYNOPSIS
    per_sample_gene_profile [options] gene_assembly sample_bam

DESCRIPTION
    It is to compute the number of read hits and correct by the 
    copy number for 16S gene sequences.
'''

posthelp='''
OUTPUT
    A tab-delimited file, of which column is gene and row is sample.

EXAMPLES
    -
 
SEE ALSO
    -

AUTHORS
    Feng Zeng
'''

def extract_bam_reads(bam, verbose=False):
    '''
    It is to extract short reads from the input bam file.
    '''
    if verbose:
        logging.info('extract short reads from ' + os.path.basename(bam))

    prefix = os.path.join(os.getcwd(),os.path.splitext(os.path.basename(bam))[0])
    PE1,PE2,SIN = extract_reads.extract_reads(bamlist=[bam],prefix=prefix)
    # merge
    os.system('cat {} {} {} | seqtk seq -A > {}'.format(PE1,PE2,SIN,prefix+'.fasta'))
    # remove temporary file
    cmd = ['rm','-f',PE1,PE2,SIN]
    call(cmd)
    return prefix+'.fasta'


def make_blast_database(gene_seq,verbose=False):
    '''
    It is to make index database for gene sequences.
    '''
    if verbose:
        logging.info('make Blast index of '+os.path.basename(gene_seq))

    cmd = ['makeblastdb','-dbtype','nucl']
    cmd += ['-in', gene_seq]
    cmd += ['-out',os.path.join(os.getcwd(),os.path.basename(gene_seq))]
    call(cmd,stderr=fnull,stdout=fnull)
    return os.path.join(os.getcwd(),os.path.basename(gene_seq))


def run_blast(gene_db, reads, word_size=22, reward=1, penalty=-2, evalue=1e-30, cores=1, verbose=False):
    '''
    It is to blast reads against genes.
    '''
    if verbose:
        logging.info('run Blast using '+str(cores)+' cores for '+os.path.basename(gene_db)+' and '+os.path.basename(reads))

    prefix = gene_db

    cmd = ['blastn']
    cmd += ['-word_size',str(word_size)]
    cmd += ['-reward',str(reward)]
    cmd += ['-penalty',str(penalty)]
    cmd += ['-num_threads',str(cores)]
    cmd += ['-evalue',str(evalue)]
    cmd += ['-db',gene_db]
    cmd += ['-query',reads]
    cmd += ['-out',prefix+'_blast.xml']
    cmd += ['-outfmt','5']
    call(cmd)
    return prefix+'_blast.xml'


def get_read_hits(blast_results,verbose):
    '''
    It is to parse the blast results and get the hit information
    for reads.
    '''
    if verbose:
        logging.info('parse Blast results')

    read_hits = defaultdict()

    with open(blast_results,'r') as f:
        for item in NCBIXML.parse(f):
            read = item.query
            read_length = int(item.query_length)
            if read.endswith('/1') or read.endswith('.1'):
                read = read[:-2]
            elif read.endswith('/2') or read.endswith('.2'):
                read = read[:-2]

            best_evalue = 100
            best_genes = []
            for hit in item.alignments:
                for hsp in hit.hsps:
                    if float(hsp.identities) < .95*float(hsp.align_length):
                        continue
                    evalue = float(hsp.expect)
                    if evalue < best_evalue:
                        best_evalue = evalue
                        best_genes = [hit.hit_def.split()[0]]
                    elif evalue == best_evalue:
                        best_genes += [hit.hit_def.split()[0]] 
            if len(best_genes)==0:
                continue
            if read in read_hits:
                if read_hits[read]['best_score'] > best_evalue:
                    read_hits[read]['best_score'] = best_evalue
                    read_hits[read]['genes'] = list(Set(best_genes))
                elif read_hits[read]['best_score'] == best_evalue:
                    read_hits[read]['genes'] = list(Set(read_hits[read]['genes']+best_genes))
            else:
                read_hits.setdefault(read,{'best_score':best_evalue,'genes':list(Set(best_genes))})

    return read_hits


def get_gene_count(read_hits):
    gene_count = defaultdict(float)
    for read, record in read_hits.iteritems():
        Z = len(record['genes'])
        for g in record['genes']:
            if g in gene_count:
                gene_count[g] += round(1./Z,5);
            else:
                gene_count[g] = round(1./Z,5);
    return gene_count
    
def gene_read_count(blast_results, verbose):
    if verbose:
        logging.info('parse blast output')

    read_hits = get_read_hits(blast_results)
    gene_count = get_gene_count(read_hits)
    return gene_count


def load_gene_copy_number():
    taxon_copy_number = defaultdict(float)
    # parse the copy number file
    with open(rrnDB) as f:
        reader = csv.DictReader(f,delimiter='\t')
        for row in reader:
            taxon_copy_number.setdefault(row['name'],float(row['mean']))
    return taxon_copy_number
    

phylogeny_level = ['domain','phylum','class','order','family','genus']
phylogeny_level_prefix = {'domain':'0__',\
                          'phylum':'1__',\
                          'class':'2__',\
                          'order':'3__',\
                          'family':'4__',\
                          'genus':'5__'}

def parse_rdp_result(rdp_result,thresh):
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
                    taxon = taxon.strip()
                    if score>=thresh:
                        lineage += [phylogeny_level_prefix[level]+taxon]
                    taxon,level,score = '','',0
                    terms = []
            gene_lineage.setdefault(gene,lineage)
    return gene_lineage
    
def copy_number_correct(gene_seq,gene_count,rdp_classifier,taxon_copy_number,thresh):
    prefix = os.path.splitext(gene_seq)[0]
    rdp_out = prefix+'_rdp_result.txt'
    # 1 RDP classification
    cmd = ['java','-jar',rdp_classifier,'-f','fixrank','-o',rdp_out,gene_seq]
    call(cmd)
    # 2 parse RDP result
    gene_lineage = parse_rdp_result(rdp_out,thresh)
    # 3 correction
    for gene,count in gene_count.iteritems():
        lineage = sorted(gene_lineage[gene])
        lineage_copy_number = []
        for i,t in enumerate(lineage):
            t = t.split('__')[1]
            if t in taxon_copy_number:
                lineage_copy_number += [taxon_copy_number[t]]
        if len(lineage_copy_number)==0:
            copy_number = 1.0
        else:
            copy_number = lineage_copy_number[-1]
        gene_count[gene] = count/copy_number  

    call(['rm',rdp_out])

    return gene_count
    

def per_sample_gene_profile(gene_seq, reads, taxon_copy_number=None, \
                            rdp_classifier=None, thresh=0.6, \
                            word_size=22, reward=1, penalty=-2, \
                            evalue=1e-30, cores=1, verbose=False):
    # 1 make blast database
    gene_db = make_blast_database(gene_seq,verbose)

    # 2 blasting
    blast_results = run_blast(gene_db=gene_db,\
                              reads=reads,\
                              word_size=word_size,\
                              reward=reward,\
                              penalty=penalty,\
                              evalue=evalue,\
                              cores=cores,\
                              verbose=verbose)

    # 3 counting
    gene_count = gene_read_count(blast_results, verbose)

    # 4 correcting
    if taxon_copy_number:
        gene_count = copy_number_correct(gene_seq,gene_count,rdp_classifier,taxon_copy_number,thresh)

    return gene_count


def per_sample_read_hit(gene_seq, reads, \
                        word_size=22, reward=1, penalty=-2, \
                        evalue=1e-30, cores=1, verbose=False):
    # 1 make blast database
    gene_db = make_blast_database(gene_seq,verbose)

    # 2 blasting
    blast_results = run_blast(gene_db=gene_db,\
                              reads=reads,\
                              word_size=word_size,\
                              reward=reward,\
                              penalty=penalty,\
                              evalue=evalue,\
                              cores=cores,\
                              verbose=verbose)

    # 3 counting
    read_hits = get_read_hits(blast_results,verbose)

    return read_hits

random_string='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'

def per_sample_read_hit_parallel_wrapper(params):
    cpu_id,cores,gene_seq,gene_set,reads = params[:5]
    word_size,reward,penalty,evalue = params[5:9]
    verbose = params[9]

    if verbose:
        logging.info('process at cpu#'+str(cpu_id))

    wd = os.getcwd()

    # make the temporary directory
    twd = os.path.join(wd,'temp_cpu_'+str(cpu_id)+'_'+''.join(random.sample(random_string,5)))
    os.mkdir(twd)

    os.chdir(twd)

    prefix = os.path.splitext(os.path.basename(gene_seq))[0]
    # extract gene sequence
    gene_chunk_seq = os.path.join(twd,prefix+'_'+str(cpu_id)+'.fasta')
    cmd = ['samtools','faidx',gene_seq] + list(gene_set)
    with open(gene_chunk_seq,'w') as f:
        call(cmd,stdout=f)
    
    # counting
    read_hits = per_sample_read_hit(gene_seq=gene_chunk_seq,\
                                    reads=reads,\
                                    word_size=word_size,reward=reward,penalty=penalty,\
                                    evalue=evalue,cores=max([1,cpu_count()/cores]),\
                                    verbose=verbose)

    # return back
    os.chdir(wd)
    shutil.rmtree(twd)
 
    return read_hits


def main():
    global opts

    wd = os.getcwd()
    twd = os.path.join(wd,'temp_gene_count_'+''.join(random.sample(random_string,5)))

    os.mkdir(twd)
    os.chdir(twd)

    # parameter setting
    gene_seq = opts.fasta
    read_bam = opts.bam

    cores = min([opts.cores,cpu_count()])
    verbose = opts.verbose

    word_size = opts.word_size
    reward = opts.reward
    penalty = opts.penalty
    evalue = opts.e_value

    ignore_copy_correct = opts.ignore_copy_correct
    rdp_classifier = opts.rdp_classifier
    thresh = opts.thresh

    # build the fasta index
    cmd = ['samtools','faidx',gene_seq]
    call(cmd)

    # split genes into chunks
    gene_set = []
    with open(gene_seq+'.fai') as f:
        for line in f:
            gene_set += [line.split()[0]]

    gene_chunk = np.array_split(gene_set,cores)

    # extract reads
    reads = extract_bam_reads(read_bam, verbose)

    # parameters for chunks
    params = []
    for i,gene_subset in enumerate(gene_chunk):
        params.extend([(i,cores,gene_seq,gene_subset,reads,word_size,reward,penalty,evalue,verbose)])

    # create multithreads
    read_hits_chunk = []
    pool = Pool(cores)
    r = pool.map_async(per_sample_read_hit_parallel_wrapper,params,callback=read_hits_chunk.extend)
    r.wait()
    pool.close()
    pool.join()

    # merge chunk counts
    read_hits = defaultdict()
    for rh in read_hits_chunk:
        for r in rh.keys():
            if r in read_hits:
                if read_hits[r]['best_score']>rh[r]['best_score']:
                    read_hits[r]['best_score'] = rh[r]['best_score']
                    read_hits[r]['genes'] = rh[r]['genes']
                elif read_hits[r]['best_score']==rh[r]['best_score']:
                    read_hits[r]['genes'] += rh[r]['genes']
            else:
                read_hits[r] = rh[r]

    # compute gene count
    gene_count = get_gene_count(read_hits)

    # copy number correct
    if not ignore_copy_correct:
        if verbose:
            logging.info('copy number correction')
        taxon_copy_number = load_gene_copy_number()
        gene_count = copy_number_correct(gene_seq,gene_count,rdp_classifier,taxon_copy_number,thresh)

    Z = 1
    if opts.relative_abundance:
        for g,c in gene_count.iteritems():
            Z += c

    # dump to local file
    fieldnames = ['sample'] + [opts.sample]
    with open(os.path.join(wd,opts.sample+"_gene_count.tsv"),'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames,delimiter='\t')
        writer.writeheader()
        for g,c in sorted([(g,c) for g,c in gene_count.items()]):
            if opts.relative_abundance:
                c = c/Z
            else:
                c = int(round(c))
            writer.writerow({'sample':g,opts.sample:c})

    os.chdir(wd)
    shutil.rmtree(twd)
    

if __name__=='__main__':
    try:
        t_start = time.time()
        # set command-line parser
        parser = argparse.ArgumentParser(description=prehelp,epilog=posthelp,\
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument('fasta',help='gene sequences',metavar='GENE_SEQ')
        parser.add_argument('bam',help='sample reads',metavar='SAMPLE_BAM')
        parser.add_argument('sample',help='sample name',metavar='SAMPLE_NAME')
        parser.add_argument('-c','--cores',help='specify the CPU cores [1]',default=1,type=int,dest='cores')
        parser.add_argument('-w','--word-size',help='word size for Blast [22]',\
                            default=22,type=int,dest='word_size')
        parser.add_argument('-R','--reward',help='reward score for Blast [1]',default=1,type=int,dest='reward')
        parser.add_argument('-P','--penalty',help='penalty score for Blast[-2]',default=-2,type=int,dest='penalty')
        parser.add_argument('-e','--e-value',help='e-value threshold for Blast [1e-10]', \
                            default=1e-10,type=float,dest='e_value')
        parser.add_argument('-n','--ignore-copy-correct',help='turn off copy number correction [False]',\
                            dest='ignore_copy_correct',\
                            action='store_true',default=False)
        parser.add_argument('-C','--rdp-classifier',help='path to the RDP classifier',dest='rdp_classifier',\
                            default='/home/fzeng/Tool/MetagenomeTools/RDPTools/classifier.jar',type=str)
        parser.add_argument('-t','--thresh',help='specify phylogeny assignment threshold for copy number correction [0.6]',\
                            dest='thresh',default=0.6,type=float)
        parser.add_argument('-r','--rel',help='output relative abundance [False]',dest='relative_abundance',\
                            action='store_true',default=False)
        parser.add_argument('-v','--verbose',help='verbose output',dest='verbose',action='store_true',default=False)

        # parse options
        opts = parser.parse_args()
        opts.fasta = os.path.abspath(opts.fasta)
        opts.bam = os.path.abspath(opts.bam)

        opts.cores = min([opts.cores,cpu_count()])

        # set logging
        logging.basicConfig(format='[%(asctime)s] %(levelname)s : %(message)s', level=logging.INFO)

        # check rdp classifier
        if not opts.ignore_copy_correct and not os.path.exists(opts.rdp_classifier):
            print >>sys.stderr,"Need to specify the RDP classifier for copy number correction!"
            print >>sys.stderr,"Exit the program!"
            sys.exit(0)

        # run the program
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
