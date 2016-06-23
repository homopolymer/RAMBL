#!/usr/bin/env python
"""
Calculate depth and breadth of marker genes across multiple samples.
"""

from __future__ import division
import os
import sys
import time
import logging
import argparse
import tempfile
import shutil
import subprocess
import itertools
from multiprocessing import Pool
import numpy as np
import pipes


def batch_depth(params):
    global opts
    tmp_dir,i,bamlist = params
    gc = os.path.join(tmp_dir,"depth_on_local_bamlist_%d.cov"%i)
    if opts.verbose:
        logging.info('calculate depth on bamlist#%d to %s' %(i,gc))
    # run samtools depth on local bamlist
    f = open(gc,'w')
    subprocess.call(['samtools','depth']+[str(x) for x in bamlist],stdout=f)
    f.close()
    subprocess.call(['sort','-k1,1n','-k2,2n','-o',str(gc),str(gc)])
    # convert depth to bed
    bed = os.path.join(tmp_dir,"depth_on_local_bamlist_%d.bed"%i)
    if opts.verbose:
        logging.info('save depth of bamlist#%d to %s' % (i,bed))
    p = pipes.Template()
    p.append(r'awk %s{a=0;for(i=3;i<=NF;i++){a+=$i}print $1"_"$2"_"($2+1)"\t"a}%s %s' % ("'","'",str(gc)),'--')
    p.append('sort -k1,1n','--')
    f = p.open(bed,'w')
    f.close()
    # clearn gc
    subprocess.call(['rm',str(gc)])
    return bed



def nucleotide_depth(params):
    global opts
    tmp_dir,i,M,bedsplit = params
    if opts.verbose:
        logging.info('summarize nucleotide-based depth on cpu%d'%i)
 
    # extract ids
    p = pipes.Template()
    p.append('cut -f1 %s' % bedsplit,'--')
    p.append(r'sed %ss/_/\t/g%s'%("'","'"),'--')
    bed_id = os.path.join(tmp_dir,'%s.id' % bedsplit)
    f = p.open(bed_id,'w')
    f.close()

    # compute sum
    p = pipes.Template()
    p.append('cut -f2-%d %s'%(1+M,bedsplit),'--')
    p.append(r'sed %ss/[[:space:]]\+/+/g%s'%("'","'"),'--')
    p.append('bc','--')
    bed_sum = os.path.join(tmp_dir,'%s.sum' % bedsplit)
    f = p.open(bed_sum,'w')
    f.close()

    # join id and sum
    p = pipes.Template()
    p.append('paste %s %s'%(bed_id,bed_sum),'--')
    #p.append('sort -k1,1n -k2,2n','--')
    f = p.open("%s.0"%bedsplit,'w')
    f.close()
    subprocess.call(['rm','%s.id'%bedsplit,'%s.sum'%bedsplit])

    if opts.verbose:
        logging.info('run bedtools-merge on cpu%d'%i)
    # bedtools merge
    p = pipes.Template()
    p.append('sort -k1,1n -k2,2n %s.0'%bedsplit,'--')
    p.append('bedtools merge -c 4 -o mean -d 10 -i stdin','--')
    f = p.open(bedsplit,'w')
    f.close()
    return bedsplit


def main():
    global opts
    try:
        # TODO: load fofn 
        with open(opts.fofn[0],'r') as f:
            bamlist = [x.rstrip() for x in f]
        N = len(bamlist)
        M = opts.cores
        bamlist = np.array_split(bamlist,M)
        if opts.verbose:
            logging.info('%d bam files in total, and split to %d parts' % (N,M))

        # TODO: build temporary directory
        tmp_dir = tempfile.mkdtemp()
        if opts.verbose:
            logging.info('work in temporary dir %s' % tmp_dir)

        # TODO: run samtools depth on local bamlist
        params = []
        for i,local_bamlist in enumerate(bamlist):
            params.extend([(tmp_dir,i,list(local_bamlist))])

        pool = Pool(opts.cores)
        bedlist = []
        r = pool.map_async(batch_depth,params,callback=bedlist.extend)
        r.wait()
        pool.close()
        pool.join()
        
        # TODO: run bedtools on bedlist
        if opts.verbose:
            logging.info('merge the results')

        # TODO: join local results
        if opts.verbose:
            Z = 0
            for b in bedlist:
                Z += int(os.path.getsize(b))
            logging.info('join {} depth-files in total size {} Mb'.format(M,round(Z/1024./1024.,5))) 

        bed = os.path.join(tmp_dir,'all.bed')
        if len(bedlist)==1:
            subprocess.call(['mv',bedlist[0],bed])
        else:
            p = pipes.Template()
            p.append(r'join -t $%s\t%s -e 0 -a 1 -a 2 -j 1 -o 0,1.2,2.2 %s %s 2>/dev/null'%\
                     ("'","'",bedlist[0],bedlist[1]),\
                      '--')
            for i in xrange(2,len(bedlist)):
                p.append(r'join -t $%s\t%s -e 0 -a 1 -a 2 -j 1 -o 0,%s,2.2 - %s 2>/dev/null'%\
                         ("'","'",','.join(['1.%d' % (x+2) for x in xrange(i)]),bedlist[i]),\
                          '--')
            f = p.open(bed,'w')
            f.close()

        # TODO: split into parts
        wc = subprocess.Popen(['wc','-l',str(bed)],stdout=subprocess.PIPE)
        wc_res,wc_err = wc.communicate()
        LN = int(wc_res.strip().split()[0])
        if opts.verbose:
            logging.info('split %d nucleotide positions into %d parts of each with %d positions' %\
                        (LN,opts.cores,int(np.ceil(LN/float(opts.cores)))))
        subprocess.call(['split','-a','3','-d','-l',str(int(np.ceil(LN/float(opts.cores)))),str(bed),\
                        str(os.path.join(tmp_dir,'all_bed_split'))])
        bedsplit = [os.path.join(tmp_dir,'all_bed_split%03d'%x) for x in xrange(opts.cores)]

        params = []
        for i,bs in enumerate(bedsplit):
            params.extend([(tmp_dir,i,M,bs)])

        pool2 = Pool(opts.cores)
        bedsplit = []      
        r2 = pool2.map_async(nucleotide_depth,params,callback=bedsplit.extend)
        r2.wait()
        pool2.close()
        pool2.join()

        # TODO: report
        if opts.verbose:
            logging.info('merge %d parts, then sort' % opts.cores)
        p = pipes.Template()
        p.append('cat %s' % " ".join([str(x) for x in bedsplit]),'--')
        p.append('sort -k1,1n -k2,2n','--')
        f = p.open(bed,'w')
        f.close()
        proc_sum_bed = subprocess.Popen(['bedtools','merge','-c','4','-o','mean','-d','10','-i',str(bed)])
        proc_sum_bed.wait()
        if opts.verbose:
            logging.info('rm %s' % tmp_dir)
        #shutil.rmtree(tmp_dir)
        if opts.verbose:
            logging.info('complete')

    except Exception as e:
        logging.exception('rm %s' % tmp_dir)
        shutil.rmtree(tmp_dir)
        raise e    


if __name__ == "__main__":
    try:
        tstart = time.time()
        # TODO: command-line parser
        parser = argparse.ArgumentParser(description=globals()['__doc__'])
        parser.add_argument('fofn', nargs=1, help="list of bam filenames, one filename by one line")
        parser.add_argument('gi', nargs=1, help="gene/genome index file")
        parser.add_argument('-c','--cores',help='computing cores [1]',type=int,default=1,metavar='INT')
        parser.add_argument('-v', action="store_true", help="verbose output", default=False, dest="verbose")
        # TODO: parsing
        opts = parser.parse_args()
        logging.basicConfig(format='[%(asctime)s] %(levelname)s : %(message)s', level=logging.INFO)
        # TODO main routine
        main()

        tescape = round((time.time()-tstart)/60.0,5)
        if opts.verbose:
            logging.info('elapsed time is {} minutes'.format(tescape))
        sys.exit(0)
    except KeyboardInterrupt as e:
        raise e
    except SystemExit as e:
        raise e
    except Exception as e:
        logging.exception('ERROR, UNEXPECTED EXCEPTION')
        logging.exception(str(e))
        traceback.print_exc()
        os._exit(1)
