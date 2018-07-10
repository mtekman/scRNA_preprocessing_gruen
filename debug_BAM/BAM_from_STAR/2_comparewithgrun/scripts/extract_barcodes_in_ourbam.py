#!/usr/bin/env python

import sys
import os.path
import pickle
import pysam
from GeneralUtils import GeneralUtils as gu

if len(sys.argv) != 5:
    print('''
Takes a BAM file with UMIs and cell barcodes already appended to the readname (e.g. via umi_tools and STAR), the wanted barcodes, and *for a given gene of interest*, performs two tasks:

 1. Scores each read with an annotation (e.g. 'f' for flagged, 'n' for NM > 2, 'x' for nm <= NM) on whether it will be used in the count matrix.

 2. Produces a matrix of UMIs against Barcodes

 This matrix shows how many PCR duplicates exist within a cell for a specific gene, and can be collapsed later to yield the number of unique transcripts for that gene in a particular cell.  Do this for all genes and merge tables, and you will yield the actual count matrix.

The gene of interest should be given by the reference coordinates as determined by the reference.

e.g. for ENSDARG00000019692, this is chr16:20392245-20433166 in GRCz10/danRer10


%s <STAR BAM> <chrJ:123456:7891011> <wanted_barcodes> <out_file>

''' % sys.argv[0],
          file=sys.stderr)
    exit(-1)
    
# Parse
wanted_bars = gu.wantedBarcodes(sys.argv[3])
chrom, loc1, loc2 = gu.parseLoci(sys.argv[2])
samfile = pysam.AlignmentFile(sys.argv[1], "rb")

outfile = open(sys.argv[4],'w')

for read in samfile.fetch(chrom, loc1, loc2):
    NM = read.get_tag('nM') # number of mismatches to reference (<= 2)
    NH = read.get_tag('NH') # number of multimapped reads (== 1)

    # flag indicating strandedness and multimapped, sometimes indicative of NH
    # but also sometimes not
    flag = read.flag

    bad_flag = flag not in [0,16]   # both forward and reverse
    bad_NM = NM > 2
    bad_NH = NH > 1

    #print(flag, read.get_tags())
    initial = "%s%s%s" % (
        "f" if bad_flag else "",
        "n" if bad_NM else "",
        "h" if bad_NH else "",
    )
    
    print(initial + read.tostring(), file=outfile)
   
    if not(bad_flag or bad_NM or bad_NH):
        gu.updateMap(read.query_name)
        

gu.printMap(sys.argv[2], outfile.name)
