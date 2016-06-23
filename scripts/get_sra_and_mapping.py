#!/usr/bin/env python
'''
Download SRA files, extract to FASTQ files, and map to gene references.

Example
=======
get_sra_and_map.py [options] SraRunTable.txt
'''

import csv
import time
import os,sys
import logging
import argparse
import subprocess
from collections import defaultdict
from multiprocessing import Pool,cpu_count

def get_sra(args):
    global cwd,opts 
    smpl,runs = args
    #ftp_prefix='ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra'
    # new ftp prefix
    ftp_prefix='ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra'
    os.chdir(os.path.join(cwd,smpl))
    for run in runs:
        sra_data = '%s/%s/%s/%s/%s.sra'%(ftp_prefix,run[:3],run[:6],run,run)
        cmd = ['curl','-C','-','-O',sra_data]
        if opts.verbose:
            logging.info('get sra file of '+run)
        sys.stderr.write(smpl+'\t'+' '.join(cmd)+'\n')
        subprocess.call(cmd)
    os.chdir(cwd)
    return None

def extract_fastq(args):
    global cwd,opts
    smpl,runs = args
    os.chdir(os.path.join(cwd,smpl))
    for run in runs:
        cmd = ['fastq-dump','--split-3','-W','--skip-technical','%s.sra'%run]
        if opts.verbose:
            logging.info('extract fastq files of '+run)
        sys.stderr.write(smpl+'\t'+' '.join(cmd)+'\n')
        subprocess.call(cmd)
    os.chdir(cwd)
    return None

def mapping(args):
    global cwd,opts
    smpl,geneseq,geneidx = args
    # mapping
    cmd = ['bowtie2','--local','-p',str(min([opts.cores,cpu_count()])),'-x',geneidx]
    os.chdir(os.path.join(cwd,smpl))
    for fq in os.listdir('.'):
        if fq.endswith('.fastq'):
            if '_1' in fq:
                cmd.extend(['-1',fq])
            elif '_2' in fq:
                cmd.extend(['-2',fq])
            else:
                cmd.extend(['-U',fq])
    cmd.extend(['-S','%s_%s.sam'%(smpl,opts.suffix)])            
    sys.stderr.write(smpl+' '+' '.join(cmd)+'\n')
    subprocess.call(cmd)
    # convert sam to bam
    cmd = ['samtools','view','-F4','-bht','%s.fai'%geneseq,'%s_%s.sam'%(smpl,opts.suffix)]
    if opts.verbose:
        logging.info('mapping reads of '+smpl)
    sys.stderr.write(smpl+' '+' '.join(cmd)+'\n')
    with open('%s_%s.bam'%(smpl,opts.suffix),'w') as f:
        subprocess.call(cmd,stdout=f)
    # sort bam
    cmd = ['samtools','sort','%s_%s.bam'%(smpl,opts.suffix),'temp']
    sys.stderr.write(smpl+' '+' '.join(cmd)+'\n')
    subprocess.call(cmd)
    # mv temp bam to final bam
    cmd = ['mv','temp.bam','%s_%s.bam'%(smpl,opts.suffix)]
    sys.stderr.write(smpl+' '+' '.join(cmd)+'\n')
    subprocess.call(cmd)
    # indexing
    cmd = ['samtools','index','%s_%s.bam'%(smpl,opts.suffix)]
    sys.stderr.write(smpl+' '+' '.join(cmd)+'\n')
    subprocess.call(cmd)
    # remove sam
    cmd = ['rm','%s_%s.sam'%(smpl,opts.suffix)]
    sys.stderr.write(smpl+' '+' '.join(cmd)+'\n')
    subprocess.call(cmd)
    os.chdir(cwd)
    return None


if __name__=="__main__":
    t_start = time.time()
    
    geneidx = '/data/fzeng/Genome/OTU/GreenGenes/gg_13_8_otus/index/bt2/99_otus.fasta'
    geneseq = '/data/fzeng/Genome/OTU/GreenGenes/gg_13_8_otus/rep_set/99_otus.fasta'

    # set command-line parser
    parser = argparse.ArgumentParser(description=globals()['__doc__'],formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('sra_table',help='data collection',metavar="SRA")
    parser.add_argument('-w','--start-with',help='''run the script starting from a step\n
    =0    run downloading SRA files\n
    =1    run extracting FASTQ files\n
    =2    run mapping to gene references\n
    >0    start from downloading\n
    >1    start from extracting\n
    <1    only run downloading\n
    <2    only run downloading and extracting''', default='>0',type=str,dest='start_with',metavar='STR')
    parser.add_argument('-c','--cores',help='number of computing cores [1]',\
                        default=1,type=int,dest='cores',metavar='INT')
    parser.add_argument('-o','--outdir',help='directory to save files',metavar='STR')
    parser.add_argument('-s','--suffix',help='suffix of mapping results [to_gg_99_otus]',\
                        default='to_gg_99_otus',type=str,metavar='STR')
    parser.add_argument('-x','--gene-dict',help='Bowtie2 index location',default=geneidx,\
                        type=str,dest='gene_dict',metavar='STR')
    parser.add_argument('-f','--gene-fa',help='gene fasta file',default=geneseq,\
                        type=str,dest='gene_fa',metavar='STR')
    parser.add_argument('-v',dest='verbose',action='store_true',default=False,help='verbose output')
    # parse options
    opts = parser.parse_args()

    # set logging format
    logging.basicConfig(format='[%(asctime)s] %(levelname)s : %(message)s', level=logging.INFO)

    if opts.outdir:
        if os.path.exists(os.path.dirname(opts.outdir)):
            cwd = os.path.abspath(opts.outdir)
        else:
            cwd = os.path.join(os.getcwd(),opts.outdir)
    else:
        cwd = os.getcwd()

    if not os.path.exists(cwd):
        os.mkdir(cwd)

    data = defaultdict(list)
    with open('SraRunTable.txt') as f:
        reader = csv.DictReader(f,delimiter='\t')
        for row in reader:
            data.setdefault(row['Sample_Name_s'],[]).extend([row['Run_s']])

    # make dir for sample
    for smpl in data.keys():
        if not os.path.exists(os.path.join(cwd,smpl)):
            os.mkdir(os.path.join(cwd,smpl))

    if opts.start_with=='=0' or opts.start_with=='>0' or opts.start_with=='<1' or opts.start_with=='<2':
        # get sra file
        map(get_sra,[(smpl,runs) for smpl,runs in data.iteritems()])

    if opts.start_with=='=1' or opts.start_with=='>0' or opts.start_with=='>1' or opts.start_with=='<2':
        # extract_fastq 
        pool = Pool(min([opts.cores,cpu_count()]))
        R = []
        r = pool.map_async(extract_fastq,[(smpl,runs) for smpl,runs in data.iteritems()],callback=R.extend)
        r.wait()
        pool.close()
        pool.join()

    if opts.start_with=='=2' or opts.start_with=='>0' or opts.start_with=='>1':
        # mapping
        map(mapping,[(smpl,opts.gene_fa,opts.gene_dict) for smpl in data.keys()])
