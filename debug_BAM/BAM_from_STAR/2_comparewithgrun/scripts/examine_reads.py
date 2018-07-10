#!/usr/bin/env python3

import sys
import pysam

if len(sys.argv) < 3:
    print('''
SAM files filtered for a specific gene or locus, and a list of umi:bar combos to extract from them and print them out in org-mode blocks.

%s <sam_grun> <sam_mine> <[umi:bar]>
''' % sys.argv[0], file=sys.stderr)
    exit(-1)


sam_file_grun=sys.argv[1]
sam_file_mine=sys.argv[2]

umi_bar_combos=sys.argv[3:]
# umi:bar umi:bar, etc


mer = {} # umi -> bar -> G: lines, M:lines

for umi_bar in umi_bar_combos:
    umi, bar = umi_bar.split(':')

    if umi not in mer:
        mer[umi] = {}

    mer[umi][bar] = {'G':[], 'M':[]}
    

def searchSAM(file, symbol):
    with open(file, 'r') as f:
        for line in f:
            tokens = line.split()
            readname = tokens[0]
            read, bar, umi = readname.split('_')

            if umi not in mer:
                continue

            if bar not in mer[umi]:
                continue
            
            mer[umi][bar][symbol].append(line)


searchSAM(sam_file_grun, 'G')
searchSAM(sam_file_mine, 'M')

# print all
for umi_bar in umi_bar_combos:
    umi, bar = umi_bar.split(':')

    grun_reads = mer[umi][bar]['G']
    mine_reads = mer[umi][bar]['M']

    print(
        '''
#+BEGIN_SRC Grun %s:%s
%s
#+END_SRC

#+BEGIN_SRC Mine %s:%s
%s
#+END_SRC
    ''' % (
        umi, bar, '\n'.join(grun_reads),
        umi, bar, '\n'.join(mine_reads)
    ))




