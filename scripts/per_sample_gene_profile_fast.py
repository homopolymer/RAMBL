#!/usr/bin/env python

import os
import sys
import time
import glob
import logging
import argparse
from subprocess import call
import extract_reads
from sets import Set
from collections import defaultdict
import csv
from multiprocessing import cpu_count
import random
from random import choice
import numpy as np
import pandas as pd
import shutil

fnull = open(os.devnull,'w')

rrnDB=os.path.join(os.path.dirname(__file__),'rrnDB_RDP.tsv')

prehelp='''
NAME
    per_sample_gene_profile_fast

VERSION
    -

SYNOPSIS
    per_sample_gene_profile_fast [options] gene_assembly sample_bam sample_name

DESCRIPTION
    It is to compute the number of read hits and correct by the 
    copy number for 16S gene sequences. It invokes a few C++ 
    routines to speed up.
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

    execdir = os.path.dirname(__file__)
    prefix = os.path.join(os.getcwd(),os.path.splitext(os.path.basename(bam))[0])
    
    cmd = [os.path.join(execdir,"extract_mapped_reads"),bam,prefix]
    call(cmd)
    
    # merge
    PE1,PE2,SIN = prefix+".1.fa",prefix+".2.fa",prefix+".single.fa"
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


def run_blast(gene_db, reads, word_size=22, reward=1, penalty=-2, evalue=1e-30, max_num_align=30, cores=1, verbose=False):
    '''
    It is to blast reads against genes.
    '''
    if verbose:
        logging.info('run Blast using '+str(cores)+' cores')

    prefix = gene_db
    prefix = 'read_contig'
    suffix = '_blast.xml'
    
    cmd = ['blastn']
    cmd += ['-word_size',str(word_size)]
    cmd += ['-reward',str(reward)]
    cmd += ['-penalty',str(penalty)]
    cmd += ['-num_threads',str(cores)]
    cmd += ['-evalue',str(evalue)]
    cmd += ['-db',gene_db]
    cmd += ['-query',reads]
    cmd += ['-out',prefix+suffix]
    cmd += ['-outfmt','5']
    cmd += ['-max_target_seqs',str(max_num_align)]
    call(cmd)
    return prefix+suffix


def blastxml2csv(blast_results,verbose=False):
    '''
    It is to convert the blast XML file to csv file.
    '''
    if verbose:
        logging.info('convert Blast XML to CSV')

    blastXML = blast_results
    blastDB = os.path.splitext(blast_results)[0]+'.db'
    blastCSV = os.path.splitext(blast_results)[0]+'.csv'

    # convert xml to sql database
    cmd = ['bigBlastParser','--max_hit','-1','--max_hsp','-1','-o',blastDB,blastXML]
    if verbose:
        logging.info(' '.join(cmd))
    call(cmd,stdout=fnull,stderr=fnull)

    # sql command
    cmd = '''#!/usr/bin/env bash
sqlite3 {} <<!
.mode csv
.output {}
select a.query_def,b.definition,c.identity,c.align_len,c.query_from,c.query_to,c.hit_from,c.hit_to,c.evalue,a.query_len
from query as a
join hit as b on a.query_id = b.query_id
join hsp as c on b.query_id = c.query_id and b.hit_id = c.hit_id;
!'''
    sqlShell = os.path.join(os.getcwd(),'sql_script.sh')
    with open(sqlShell,'w') as f:
        f.write(cmd.format(blastDB,blastCSV)+'\n')

    call(['bash',sqlShell])

    return blastCSV


def gene_raw_abundance(blast_results,max_align_iden,evalue,verbose):
    '''
    It is to parse the blast results and get the raw abundances
    for contigs.
    '''
    blast_results = blastxml2csv(blast_results,verbose)

    if verbose:
        logging.info('compute raw abundances of contigs')

    execdir = os.path.dirname(__file__)
    # run blastout2abundance
    cmd = [os.path.join(execdir,'blastout2abundance'),'-I',str(max_align_iden),'-E',str(evalue),blast_results]
    outfile = os.path.splitext(blast_results)[0]+'.raw'
    with open(outfile,'w') as f:
        call(cmd,stdout=f)

    gene_abundances = pd.read_csv(outfile,sep='\t',header=None)
    
    return gene_abundances.set_index(0).to_dict()[1]



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
    
def copy_number_correct(gene_seq,gene_count,rdp_classifier,taxon_copy_number,thresh,verbose=False):
    if verbose:
        logging.info('normalize contig abundance by copy number')

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
                            evalue=1e-30, max_num_align=30, max_align_iden=95,\
                            cores=1, verbose=False):
    # 1 make blast database
    gene_db = make_blast_database(gene_seq,verbose)

    # 2 blasting
    blast_results = run_blast(gene_db=gene_db,\
                              reads=reads,\
                              word_size=word_size,\
                              reward=reward,\
                              penalty=penalty,\
                              evalue=evalue,\
                              max_num_align=max_num_align,\
                              cores=cores,\
                              verbose=verbose)

    # 3 counting
    gene_count = gene_raw_abundance(blast_results, max_align_iden, evalue, verbose)

    # 4 correcting
    if taxon_copy_number is not None:
        gene_count = copy_number_correct(gene_seq,gene_count,rdp_classifier,taxon_copy_number,thresh,verbose)

    return gene_count


def main():
    global opts

    wd = os.getcwd()
    twd = os.path.join(wd,'tempdir_gene_count_'+os.path.splitext(os.path.basename(opts.bam))[0])

    call(['mkdir','-p',twd])
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
    max_num_align = opts.max_num_align
    max_align_iden = opts.max_align_iden

    copy_correct = not opts.ignore_copy_correct
    rdp_classifier = opts.rdp_classifier
    thresh = opts.thresh

    # build the fasta index
    cmd = ['samtools','faidx',gene_seq]
    call(cmd)

    # extract reads
    reads = extract_bam_reads(read_bam, verbose)

    # compute gene profile
    taxon_copy_number = None
    if copy_correct:
        taxon_copy_number = load_gene_copy_number()
    gene_count = per_sample_gene_profile(gene_seq,reads,taxon_copy_number,rdp_classifier,thresh,
                                         word_size,reward,penalty,evalue,max_num_align,max_align_iden,
                                         cores,verbose)

    # convert to DataFrame
    gene_count = pd.DataFrame.from_dict({opts.sample:gene_count})
    gene_count.index.name = 'sample'

    # relative abundance
    if opts.relative_abundance:
        gene_count = gene_count/gene_count.sum()
        gene_count[opts.sample] = gene_count[opts.sample].apply(lambda x:round(x,9))
    else:
        gene_count[opts.sample] = gene_count[opts.sample].apply(lambda x:round(x,3))
    
    # write to disk file
    gene_count.to_csv(os.path.join(wd,opts.sample+"_gene_count.tsv"),sep='\t')

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
        parser.add_argument('-A','--max_num_alignments',help='maximum number of alignments for a query [30]',\
                            default=30,type=int,dest='max_num_align')
        parser.add_argument('-I','--max_identity',help='maximum alignment identity [95]',\
                            default=95,type=int,dest='max_align_iden')
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
