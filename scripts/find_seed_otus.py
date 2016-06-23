#!/usr/bin/env python
"""
Find seed OTUs based on the phylogeny tree and abundance profile

Change Log
==========
Sep 25, 2015    Feng Zeng    Add coverage criterion
Jan 10, 2016    Feng Zeng    Add clade-specific coverage
Jan 11, 2016    Feng Zeng    Add tree pruning

"""

import os
import sys
import argparse
import logging
import textwrap
import traceback
import time
import csv
import itertools
from sets import Set
from collections import defaultdict
from numpy import array,zeros
from itertools import product,combinations

from ete2 import Tree

import numpy as np
from scipy.cluster import hierarchy


from skbio.io import read as fasta_read
from skbio import DNA


gene_set_dt = np.dtype([('gene','S1024'),('abundance','f8')])



def load_gene_abundance(filename):
    '''load gene abundances from the given file

    Parameters
    ----------
    filename : gene abundance profile file

    Returns
    -------
    dict
        a dictionary object storing gene abundance

    '''
    gene_abun = defaultdict(float)
    with open(filename) as f:
        r = csv.reader(f,delimiter='\t')
        for row in r:
            gene_abun[row[0]] = float(row[3])
    return gene_abun

def load_gene_coverage(filename):
    '''load gene coverage from the given file

    Parameters
    ----------
    filename : gene abundance profile file

    Returns
    -------
    dict
        a dictionary object storing gene coverage

    '''
    gene_cover = defaultdict(float)
    with open(filename) as f:
        r = csv.reader(f,delimiter='\t')
        for row in r:
            gene_cover[row[0]] = float(row[4])
    return gene_cover


def load_gene_mask(coverage_file,gene_file):
    '''construct mask for the segment covered

    Parameters
    ----------
    coverage_file : gene coverage bed file
    gene_file     : gene faidx file

    Returns
    -------
    dict
        a dictionary object storing the mask for the segment covered
    '''
    gene_mask = defaultdict(array)  
    # parse gene index file
    with open(gene_file) as f:
        for line in f:
            fields = line.split()
            gene_mask.setdefault(fields[0],zeros(int(fields[1])))

    # parse coverage bed file
    with open(coverage_file) as f:
        for line in f:
            fields = line.split()
            p0 = int(fields[1])-1
            p1 = int(fields[2])
            gene_mask[fields[0]][p0:p1] = 1
    return gene_mask

def load_gene_size(gene_file):
    '''load gene size
    '''
    gene_size = defaultdict(int)
    with open(gene_file) as f:
        for line in f:
            fields = line.split()
            gene_size.setdefault(fields[0],int(fields[1]))
    return gene_size


def load_gene_alignment(filename):
    '''load multiple sequence alignment from the given file
    '''
    aligns = {}
    for a in fasta_read(filename,format='fasta'):
        aligns[a.metadata['id']] = DNA(a)
    return aligns


def load_phylogeny_tree(filename):
    '''load phylogeny tree from the given newick file

    Paramters
    ---------
    filename : phylogene tree file

    Returns
    -------
    tree
        a tree object

    '''
    tree = Tree(filename)
    
    # add a numeric id to tree node
    count = 0
    for node in tree.traverse('postorder'):
        node.add_feature('id',count)
        count += 1

    return tree


def remove_null_nodes(tree):
    '''delete nodes of zero abundance
    
    Parameters
    ----------
    tree : input tree object

    Returns
    -------
    tree
        a tree object
    '''
    # add remove state
    for t in tree.iter_leaves():
        if t.gene_set['abundance']>0:
            t.add_feature('prunable',False)
        else:
            t.add_feature('prunable',True)

    # postorder traversal
    for t in tree.traverse('postorder'):
        if t.is_leaf():
            continue
        child_state = [c.prunable for c in t.get_children()]
        prunable = True
        for c in t.get_children():
            prunable &= c.prunable 
        t.add_feature('prunable',prunable)

    # prunning tree
    for t in tree.traverse('postorder'):
        if t.prunable:
            t.delete()

    return tree

from itertools import izip
from random import choice as random_choice
def gene_distance(A,B):
    '''compute sequence distance between two genes A and B
    '''
    X,Y = '','' # new sequence removing common gaps
    for a,b in izip(A.values,B.values):
        if (a in A.gap_chars) and (b in B.gap_chars):
            continue
        if a in A.degenerate_chars:
            X += random_choice(list(A.degenerate_map[a]))
        else:
            X += a
        if b in B.degenerate_chars:
            Y += random_choice(list(B.degenerate_map[b]))
        else:
            Y += b
    newA = DNA(X,metadata={})
    newB = DNA(Y,metadata={})
    return newA.distance(newB)

        
def gene_tree_cluster(tree,gene_align,dissim_thres):
    '''
    cluster tree nodes based on sequence similarity
    '''
    global opts

    if tree.is_leaf():
        tree.add_feature('centroid',tree.gene_set)
        tree.add_feature('merged',1)
    else: 
        # get the children
        children = tree.get_children()

        # add @Jan 11
        # get child state
        child_merged = [child.merged for child in children]
        if np.sum(np.asarray(child_merged)==0) > 0:
            tree.add_feature('merged',0)
            return None

        # get child cluster size
        child_cluster_size = [len(child.centroid) for child in children]

        # skip if a child has multiple clusters
        if np.sum(np.asarray(child_cluster_size)>1) > 0:
            centroid = np.hstack(tuple(child.centroid for child in children))
            tree.add_feature('centroid',centroid)
            return None
        
        # create a zero distance matrix
        cluster_size = np.sum(child_cluster_size)
        dist = np.zeros((cluster_size,cluster_size))

        # fill out the up diagonal matrix
        for c1,c2 in itertools.combinations(np.arange(len(children)),2):
            d1 = int(np.sum(child_cluster_size[:c1]))  # position shift
            d2 = int(np.sum(child_cluster_size[:c2]))  # position shift
            k1 = 0
            k2 = 0
            # centroid in cluster k1 of node c1
            g1 = children[c1].centroid[k1]['gene']
            # centroid in cluster k2 of node c2
            g2 = children[c2].centroid[k2]['gene']
            # calculate the distance between cluster k1 and cluster k2
            dist_c1_c2 = tree.get_distance(str(g1),str(g2))
            #dist_c1_c2 = gene_distance(gene_align[str(g1)],gene_align[str(g2)])
            # save to distance matrix
            dist[d1+k1][d2+k2] = dist_c1_c2

        # preform hierarchical clustering on the distance matrix
        Z = hierarchy.linkage(dist)
        C = hierarchy.fcluster(Z,t=dissim_thres,criterion='distance')-1  # make sure starting from 0

        # update tree data
        centroid = np.array([(None,None)]*len(np.unique(C)),dtype=gene_set_dt)

        for c in xrange(len(children)):
            d = int(np.sum(child_cluster_size[:c]))  # position shift
            for i in xrange(child_cluster_size[c]):
                k = C[d+i]
                if np.isnan(centroid[k]['abundance']):
                    centroid[k] = children[c].centroid[i]
                else:
                    if centroid[k]['abundance'] < children[c].centroid[i]['abundance']:
                        centroid[k]['gene'] = children[c].centroid[i]['gene']
                    centroid[k]['abundance'] += children[c].centroid[i]['abundance']

        tree.add_feature('centroid',centroid)

        # merge state
        if len(np.unique(C))==1:
            tree.add_feature('merged',1)
        else:
            tree.add_feature('merged',0)


def load_gene_taxonomy(filename):
    '''load gene taxonomy annotation
    
    Parameters
    ----------
    filename : taxonomy annotation file
    
    Returns
    -------
    Dict
        a dictonary object storing gene taxonomy annotation
    '''
    gene_tax = defaultdict()
    with open(filename,'r') as f:
        for line in f:
            gene,tax = line.rstrip().split('\t')
            gene_tax[gene] = tax
    return gene_tax


def find_seed_otus():
    '''find seed otus based on phylogeny and gene abundance

    Parameters
    ----------

    Returns
    -------

    '''
    global opts

    # TODO load phylogeny tree
    if opts.verbose:
        logging.info('load phylogeny tree file')
    t0 = time.time()
    gene_tree = load_phylogeny_tree(os.path.abspath(opts.tree[0]))
    if opts.verbose:
        logging.info('elapsed time on loading phylogeny tree is {} minutes'.format(round((time.time()-t0)/60.,5)))

    # TODO load gene abundance
    if opts.verbose:
        logging.info('load gene abundance file')
    t0 = time.time()
    gene_abun = load_gene_abundance(os.path.abspath(opts.abun[0]))
    if opts.verbose:
        logging.info('elapsed time on loading gene abundances is {} minutes'.format(round((time.time()-t0)/60.,5)))

    # TODO load gene coverage
    gene_cover = load_gene_coverage(os.path.abspath(opts.abun[0]))

    # TODO load gene portion mask
    gene_mask = load_gene_mask(os.path.abspath(opts.mask[0]),os.path.abspath(opts.seq[0]))

    # TODO load gene size
    gene_size = load_gene_size(os.path.abspath(opts.seq[0]))

    # TODO load gene alignment
    gene_align = load_gene_alignment(os.path.abspath(opts.align[0]))
   
    # TODO load taxonomy annotation
    if opts.taxonomy is not None:
        if opts.verbose:
            logging.info('load gene taxonomy annotation')
        t0 = time.time()
        gene_tax = load_gene_taxonomy(os.path.abspath(opts.taxonomy))
        if opts.verbose:
            logging.info('elapsed time on loading gene taxonomy annotation is {} minutes'.format(\
                         round((time.time()-t0)/60.,5)))

    # TODO add new attributes to tree nodes
    for node in gene_tree.iter_leaves():
        if node.name in gene_abun:
            node.add_feature('gene_set',np.array([(node.name,gene_abun[node.name])],dtype=gene_set_dt))
        else:
            node.add_feature('gene_set',np.array([(node.name,0)],dtype=gene_set_dt))

    # TODO remove zero-abundance gene from tree
    if opts.verbose:
        logging.info('remove zero-abundance genes from tree')
    t0 = time.time()
    gene_tree = remove_null_nodes(gene_tree)
    if opts.verbose:
        logging.info('elapsed time on tree pruning is {} minutes'.format(round((time.time()-t0)/60.,5)))

    # TODO propogate gene cluster from leaves to root
    if opts.verbose:
        logging.info('propogate gene cluster from bottom to top based on similarity')
    t0 = time.time()
    for t in gene_tree.traverse('postorder'):
        gene_tree_cluster(t,gene_align,1.-opts.sim_thres[0])
    if opts.verbose:
        logging.info('elapsed time on propogating gene clustering is {} minutes'.format(round((time.time()-t0)/60.,5)))
   
    # TODO compute total abundance
    total_abun = 0
    for node in gene_tree.iter_leaves():
        total_abun += gene_abun[node.name]
    abun_thres = max([opts.depth_thres[0],0.0001*total_abun])
    if opts.depth_ratio is not None:
        abun_thres = opts.depth_ratio*total_abun

    # TODO report seed OTUs
    if opts.verbose:
        logging.info('find gene cluster with depth threshold {} and coverage threshold {}'.format(\
                     round(abun_thres,5),opts.gene_cover[0]))

    # save results
    seed_genes = defaultdict(dict)

    visited = Set()
    for t in gene_tree.traverse('preorder'):
        if t.id not in visited:
            visited.add(t.id)
            if t.merged == 1: # change @Jan 11
            #if len(t.centroid) == 1:
                if t.centroid[0]['abundance'] >= abun_thres:
                    gs = list()
                    gm = zeros(1)
                    if t.name in gene_cover:
                        gs.extend([(gene_cover[t.name],t.name)])
                    if t.name in gene_mask:         # Jan 10
                        gm = gene_mask[t.name]      # Jan 10
                    for n in t.get_descendants():
                        if n.name in gene_cover:
                            gs.extend([(gene_cover[n.name],n.name)])
                        if n.name in gene_mask:     # Jan 10
                            gm_l = max([gm.shape[0],gene_mask[n.name].shape[0]])  # Jan 10
                            tgm = np.resize(gm,gm_l)   # Jan 10
                            ngm = np.resize(gene_mask[n.name],\
                                            gm_l)      # Jan 10
                            gm = tgm + ngm             # Jan 10
                        visited.add(n.id)
                    gs = sorted(gs,reverse=True)
                    cover_frac = sum(gm>0)/(len(gm)+0.)        # Jan 10
                    if cover_frac >= opts.gene_cover[0]:       # Jan 10
                        g = gs[0][1]
                        tax = ''
                        if opts.taxonomy is not None:
                            tax = gene_tax[t.centroid[0]['gene']]
                        seed_genes.setdefault(g,{'gene':g,\
                                                 'clade_abundance':t.centroid[0]['abundance'],\
                                                 'clade_coverage':cover_frac,\
                                                 'seed_abundance':gene_abun[g],\
                                                 'seed_coverage':gene_cover[g],\
                                                 'taxonomy':tax,\
                                                 'count':len(t.get_leaves())})

    for g in seed_genes.keys():
        print('%(gene)s\t%(clade_abundance)f\t%(clade_coverage)f\t%(seed_abundance)f\t%(seed_coverage)f\t%(count)d\t%(taxonomy)s' % seed_genes[g])

    #seed_gene_names = sorted(seed_genes.keys())
    #seed_gene_dist = [gene_distance(gene_align[i],gene_align[j]) for i,j in combinations(seed_gene_names,2)]
    #seed_gene_z = hierarchy.linkage(seed_gene_dist,method='average')
    #seed_gene_C = hierarchy.fcluster(seed_gene_z,t=min([0.07,1-opts.sim_thres[0]]),criterion='distance')
    #if True:
    #    for c in np.unique(seed_gene_C):
    #        c_seed_gene = []
    #        for i in np.where(seed_gene_C==c)[0]:
    #            c_seed_gene += [(gene_cover[seed_gene_names[i]],seed_gene_names[i])]
    #        c_seed_gene = sorted(c_seed_gene,reverse=True)
    #        g = c_seed_gene[0][1]
    #        print('%(gene)s\t%(abundance)f\t%(depth)f\t%(ratio)f\t%(taxonomy)s' % seed_genes[g])


if __name__ == "__main__":
    try:
        # TODO set logging
        logging.basicConfig(format="[%(asctime)s] %(levelname)s : %(message)s", level=logging.INFO)

        # TODO set command-line parser
        parser = argparse.ArgumentParser(description=globals()['__doc__'],epilog=textwrap.dedent('''\
                                         Examples\n\
                                         --------\n\
                                         1. find_seed_otus 99_otus_unannotated.tree gg_99_otus_abun.txt\n\
                                         \n'''),formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument('tree',help='phylogeny tree',nargs=1,metavar="TREE")
        parser.add_argument('abun',help='abundance profile',nargs=1,metavar="ABUNDANCE")
        parser.add_argument('mask',help='gene coverage mask',nargs=1,metavar='MASK')
        parser.add_argument('seq',help='gene index file',nargs=1,metavar='GENE_INDEX')
        parser.add_argument('align',help='gene alignment file',nargs=1,metavar='GENE_ALIGN')
        parser.add_argument('-s',help='sequence similarity threshold (default: 0.9)',nargs=1,\
                            default=[0.9],type=float,dest='sim_thres',metavar="SIMILARITY")
        parser.add_argument('-d',help='depth cutting threshold (default:10)',nargs=1,\
                            default=[10],type=float,dest='depth_thres',metavar='DEPTH')
        parser.add_argument('-c',help='gene coverage threshold (default:0.6)',nargs=1,\
                            default=[0.6],type=float,dest='gene_cover',metavar='COVERAGE')
        parser.add_argument('-r',help='use adaptive depth threshold (0<=r<=1 of total depth)',\
                            type=float,dest='depth_ratio',metavar='RATIO')
        parser.add_argument('-T',help='gene taxonomy annotation',default=None,type=str,\
                            dest='taxonomy',metavar="TAXONOMY")
        parser.add_argument('-v',help='verbose optput',action="store_true",default=False,dest='verbose')
        opts = parser.parse_args()
        
        # TODO main code
        t_start = time.time()
        find_seed_otus()   

        # TODO complete
        if opts.verbose:
            logging.info('total elapsed time is {} minutes'.format(round((time.time()-t_start)/60.,5)))
        sys.exit(0)
    except KeyboardInterrupt,e:
        raise e
    except SystemExit,e:
        raise e
    except Exception,e:
        logging.exception('ERROR, UNEXCEPTED EXCEPTION')
        logging.exception(e)
        traceback.print_exc()
        os._exit(1)
