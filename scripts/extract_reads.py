#!/usr/bin/env python
'''Extract out reads from a BAM file'''

import os
import sys
import time
import logging
import argparse
import textwrap
import traceback
import subprocess
import itertools
import numpy as np
from multiprocessing import Pool
from collections import defaultdict
import glob


write_devnull = open(os.devnull,'w')

rc = {'A':'T','a':'t',\
      'C':'G','c':'g',\
      'G':'C','g':'c',\
      'T':'A','t':'a',\
      'N':'N','n':'n'}

def read_is_paired(flag):
    return ((flag&int('0x1',0))>0)

def read_is_unmapped(flag):
    return ((flag&int('0x4',0))>0)

def paired_is_unmapped(flag):
    return ((flag&int('0x8',0))>0)

def read_is_reverse_complement(flag):
    return ((flag&int('0x10',0))>0)

def read_is_paired_first(flag):
    return ((flag&int('0x40',0))>0)

def reverse_complement(seq):
    rseq = [rc[b] for b in seq[::-1]]
    return ''.join(rseq)
    
def view_bam(bam, roi=None):
    '''view a bam file with/without roi'''
    logging.debug('view %s' % bam)
    command = ['samtools','view','-F4']
    
    command.append(bam)
    
    if roi is not None and len(roi)>0:
        if isinstance(roi,(list,tuple)):
            command.extend(roi)
        else:
            command.append(roi)

    return subprocess.Popen(command,stdout=subprocess.PIPE,stderr=write_devnull).stdout

def write_read_to_file(handle,name,seq,qual):
    handle.write('@%s\n'%name)
    handle.write('%s\n'%seq)
    handle.write('+\n')
    handle.write('%s\n'%qual)


def extract_reads_of_a_bamlist(bamlist,roi=None,hs=None,h1=None,h2=None):
    '''extract sequence reads from a list of bam files'''
    logging.info('subprocess: extract reads from %d bam files' % len(bamlist))
    handles = itertools.imap(lambda b,r:view_bam(b,r),\
                             bamlist,itertools.repeat(roi))
    read_pool = defaultdict()
    # output reads
    for read in itertools.chain.from_iterable(handles):
        item = read.rstrip().split()       
        paired  = read_is_paired(int(item[1]))
        mapped  = not read_is_unmapped(int(item[1]))
        paired_mapped = not paired_is_unmapped(int(item[1]))
        if read_is_reverse_complement(int(item[1])):
            seq = reverse_complement(item[9])
            qual = ''.join([x for x in item[10][::-1]])
        else:
            seq = item[9]
            qual = item[10]
        if paired:
            read_pool.setdefault(item[0],{}).update({'type':'paired'})
            if read_is_paired_first(int(item[1])):
                read_pool.setdefault(item[0],{}).update({'1':('%s/1'%item[0],seq,qual)})
            else:
                read_pool.setdefault(item[0],{}).update({'2':('%s/2'%item[0],seq,qual)})
        else:
            read_pool.setdefault(item[0],{}).update({'type':'single'})
            read_pool.setdefault(item[0],{}).update({'data':('%s'%item[0],seq,qual)})
    for rn in sorted(read_pool.keys()):
        if read_pool[rn]['type'] == 'single':
            a,b,c = read_pool[rn]['data']
            write_read_to_file(hs,a,b,c)
        else:
            if '1' in read_pool[rn] and '2' in read_pool[rn]:
                a,b,c = read_pool[rn]['1']
                write_read_to_file(h1,a,b,c)
                a,b,c = read_pool[rn]['2']
                write_read_to_file(h2,a,b,c)
            elif '1' in read_pool[rn]:
                a,b,c = read_pool[rn]['1']
                write_read_to_file(hs,a,b,c)
            elif '2' in read_pool[rn]:
                a,b,c = read_pool[rn]['2']
                write_read_to_file(hs,a,b,c)
            

def sub_extract_reads(bamlist=list(),roilist=None,core=None,prefix=None):
    hs = open('%s.%d.single.fastq'%(prefix,core),'w')
    h1 = open('%s.%d.pe1.fastq'%(prefix,core),'w')
    h2 = open('%s.%d.pe2.fastq'%(prefix,core),'w')
    extract_reads_of_a_bamlist(bamlist,roilist,hs,h1,h2)
    hs.close()
    h1.close()
    h2.close()
    return 0
    
def sub_extract_reads_wrapper(args):
    return sub_extract_reads(*args)

def extract_reads(bamlist=list(),roilist=None,cores=1,prefix=None):
    logging.info('dump reads from %d bam files using %d cores'%(len(bamlist),min(len(bamlist),cores)))

    t0 = time.time()
    
    N = min(len(bamlist),cores)
    R = []
    pool = Pool(N)
    r = pool.map_async(sub_extract_reads_wrapper,\
                       [([bamlist[i]],roilist,i,prefix) for i in xrange(len(bamlist))],callback=R.extend)
    r.wait()
    pool.close()
    pool.join()
    
    # merge fastq files
    PE1 = sorted([f for f in glob.glob(prefix+'*.pe1.fastq')])
    cmd = ['cat'] + PE1
    with open(prefix+'_1.fastq','w') as f:
        subprocess.call(cmd,stdout=f)
    PE2 = sorted([f for f in glob.glob(prefix+'*.pe2.fastq')])
    cmd = ['cat'] + PE2
    with open(prefix+'_2.fastq','w') as f:
        subprocess.call(cmd,stdout=f)
    SIN = sorted([f for f in glob.glob(prefix+'*.single.fastq')])
    cmd = ['cat'] + SIN
    with open(prefix+'.fastq','w') as f:
        subprocess.call(cmd,stdout=f)
    
    # remove temporary files
    for f in PE1+PE2+SIN:
        subprocess.call(['rm','-f',f])

    logging.info('elapsed time on dumping reads from %d bam files is %f minutes' % (len(bamlist),(time.time()-t0)/60.))

    return (prefix+'_1.fastq',prefix+'_2.fastq',prefix+'.fastq')
    

if __name__ == "__main__":
    try:
        # set command-line parser
        parser = argparse.ArgumentParser(description=globals()['__doc__'],epilog=textwrap.dedent('''\
                                         Examples\n\
                                         --------\n\
                                         1. extract_reads\n\
                                         \n'''),formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument('bam',help='input BAM file',type=str,metavar='BAM')
        parser.add_argument('-p','--prefix',help='prefix of output read files, default: the same as input BAM file',\
                            type=str,metavar='STR')
        parser.add_argument('-v','--verbose',help='verbose output',\
                            default=False, action='store_true',dest='verbose')
        opts = parser.parse_args()

        # set logging
        if opts.verbose:
            logging.basicConfig(format="[%(asctime)s] : %(levelname)s : %(message)s", level=logging.INFO)

        # set timer
        t0 = time.time()

        # main code
        if not os.path.exists(os.path.abspath(opts.bam)):
            print >>sys.stderr,'The input BAM file is not existed.'
            sys.exit(0)
        bamlist = [os.path.abspath(opts.bam)]
           
        if opts.prefix is None:
            prefix = os.path.splitext(os.path.basename(opts.bam))[0]
        else:
            prefix = opts.prefix

        extract_reads(bamlist, None, 1, prefix)

        # complete
        logging.info('total elapsed time is %f minutes' % ((time.time()-t0)/60.))
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

    
