#!/usr/bin/env bash

mkdir -p GreenGenes
cd GreenGenes

# 16S raw sequences
curl -C - -O ftp://greengenes.microbio.me/greengenes_release/gg_13_8_otus/rep_set/99_otus.fasta

# 16S alignments
curl -C - -o 99_otus_aligned.fasta ftp://greengenes.microbio.me/greengenes_release/gg_13_8_otus/rep_set_aligned/99_otus.fasta

# 16S taxonomy
curl -C - -O ftp://greengenes.microbio.me/greengenes_release/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt

# 16S phylogenetic tree
curl -C - -O ftp://greengenes.microbio.me/greengenes_release/gg_13_8_otus/trees/99_otus_unannotated.tree 
