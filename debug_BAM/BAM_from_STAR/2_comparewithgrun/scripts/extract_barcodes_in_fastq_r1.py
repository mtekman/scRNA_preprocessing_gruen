#!/usr/bin/env python

import sys
import os.path
import pickle
from GeneralUtils import GeneralUtils as gu


if len(sys.argv)!=6:
    print('''
Using Grün's SAM file, the original FASTQ forward (barcodes+umis) file, and the wanted barcodes, *for a specific gene of interest*, performs two tasks:

 1. Appends UMIs and Cell barcodes to the readname in the SAM files, using the original barcode information from the FASTQ forward read, and also adds a score to the read depending on whether it is used in the counting process (e.g. 'f' for flagged, 'n' for NM > 2, 'x' for nm <= NM)

 2. Produces a matrix of UMIs against Barcodes

   This matrix shows how many PCR duplicates exist within a cell for a specific gene, and can be collapsed later to yield the number of unique transcripts for that gene in a particular cell.  Do this for all genes and merge tables, and you will yield the actual count matrix.

%s <Grün SAM> <FASTQ_R1> <gene of interest> <wanted_barcodes> <out_file>
''' % sys.argv[0],
          file=sys.stderr)
    exit(-1)

sam = sys.argv[1]
fastq_R1 = sys.argv[2]
gene_of_interest = sys.argv[3]
wanted_barcodes = sys.argv[4]
outfilename = sys.argv[5]

outfile = open(outfilename,'w')

wanted_bars = gu.wantedBarcodes(wanted_barcodes)

map = {}

pick = "sam."+gene_of_interest+".pickle"
if os.path.exists(pick):
    print("Loading from pickle", pick, file=sys.stderr)
    with open(pick,'rb') as pb:
        map = pickle.load(pb)

else:
    with open(sam, 'r') as f:

        for line in f:

            line = line.splitlines()[0].strip()
            tokens = line.split('\t')

            read_name = tokens[0].strip()
            gene = tokens[2].strip()

            if gene != gene_of_interest:
                continue

            if read_name in map:
                print("Duplicate", read_name, file=sys.stderr)
            else:
                map[read_name] = tokens
                print(len(map), end="\r", file=sys.stderr)

        print("Saving to pickle", pick, file=sys.stderr)
        with open(pick,'wb') as pb:
            pickle.dump(map,pb)



lc = 0
with open(fastq_R1, 'r') as f:
    for line in f:
        lc += 4

        if (lc % 400000==0):
            print(lc, end="\r", file=sys.stderr)

        name = line[1:].split()[0]
        
        sequence = f.readline()
        jj = f.readline()
        sc = f.readline()

        if name not in map:
            continue

        # Wanted read header
        # break apart umi and barcode
        umi = sequence[:6]
        bar = sequence[6:12]

        tokens = map[name]

        # Check flag
        flag = int(tokens[1])
        bad_flag = flag != 0

        # Check NM <= 2
        NM = int([t for t in tokens[10:] if t.startswith("NM:i:")][0].split('NM:i:')[-1])
        bad_NM = NM > 2

        # Check XA nm <= NM
        bad_XA = False
        if not bad_NM:
            find_XA = [t for t in tokens[10:] if t.startswith("XA:Z")]
            if len(find_XA) > 0:
                XAs = find_XA[0].split('XA:Z:')[-1].split(';')
                for xa in XAs:
                    xticks = xa.split(',')
                    if len(xticks) == 4:
                        junk, junk, junk, nm = xticks
                        if int(nm) <= NM:
                            bad_XA = True
                            break

        initial = "%s%s%s" % (
            "f" if bad_flag else "",
            "n" if bad_NM else "",
            "x" if bad_XA else ""
        )

        # Edit existing read_name
        readname = tokens[0]
        tokens[0] = initial + tokens[0] + "_" + bar + "_" + umi
        print("attached ", bar, umi, " to ", name, end="\r", file=sys.stderr)
        print('\t'.join(tokens), file=outfile)

        if not(bad_flag or bad_NM or bad_XA):
            gu.updateMap(readname, umi, bar)


gu.printMap(gene_of_interest, outfilename)

