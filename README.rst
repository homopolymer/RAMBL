****
RAMBL: Large-scale 16S gene assembly using metagenomic shotgun sequences
****

Created by Feng Zeng Summer 2015.

Copyright 2015 Feng Zeng. All rights reserved.

========
Overview
========

RAMBL is a tool for the assembly of full-length 16S genes in metagenomic shotgun data. RAMBL devices a taxonomic tree search to partition data into subgroups that comprise of short reads yielded from similar strain sequences. Within each subgroup, RAMBL deploys a progressive Dirichlet Process mixture clustering method to reconstruct full-length 16S gene sequences. 

==========
Dependency
==========


Following tools are all assumed installed and could be found in system environment::

* Python >= 2.7.10
* GCC >= 4.9.1
* Bowtie2 >= 2.2.4
* Samtools >= 1.2
* Bedtools >= 2.24.0
* Sickle >= 1.33
* Seqtk >= 1.0

Following Python packages are all assumed installed::

* ETE2 Toolkit >= 2.3.1
* Scipy >= 0.15.1
* Numpy >= 1.9.2
* Scikit-bio == 0.4.0
* Biopython >= 1.65


============
Installation
============

* Set up environment variable::
  
    $ export CXX=<PATH TO G++4.9>

* Clone the repository::

    $ git clone https://github.com/homopolymer/RAMBL.git

* Build and install to $PWD/bin::

    $ python setup.py install --prefix=$PWD

============
16S Database 
============

RAMBL relies on GreenGenes, which is available at ftp://greengenes.microbio.me/greengenes_release/gg_13_8_otus/.  The rrnDB file that records gene copy number for GreenGenes is contained in the RAMBL package.

To download GreenGenes::
    
    $ bash download_greengenes.sh


=====
Usage
=====

::

    $ python rambl.py [-c cores] [-v] [-p output_prefix] data_info.txt


=============
Data Metainfo
=============

In the above, the argument 'data_info.txt' is a metadata file.  It contains the following contents:

* Sequencing data
    1) BamFiles=<FILE>, a file listing the mapping files of data against 16S rRNA reference, one BAM file per line. ::

        $ ls -A1 *_to_gg_99_otus.bam | xargs realpath > bams.fofn

* Gene data
    1) GeneSeq=<FILE>, a FASTA file storing the reference sequences of 16S rRNA genes.
    2) GeneIndex=<FILE>, the FASTA index file.
    3) GeneTree=<FILE>, a Newick file recording the phylgenetic tree of 16S rRNA genes.
    4) GeneTax=<FILE>, a taxonomic annotation of 16S rRNA genes.
    5) GeneAlign=<FILE>, a FASTA file of 16S sequence alignments.


====
Test
====

The repository `RAMBL_Test <http://github.com/homopolymer/RAMBL_Test/>`_ contains data and script to test whether RAMBL installs and works correctly. Download test data and run the following codes for test::

    $ cd test
    $ python make_datainfo.py
    $ python <PATH TO rambl.py> -c 20 -v data_info.txt 2>&1 | tee log.txt
    $ less log.txt


================
Development Team
================

* Feng Zeng, Xiamen University
* Zicheng Wang, Tsinghua University
* Ying Wang, Xiamen University
* Ting Chen, Tsinghua University

