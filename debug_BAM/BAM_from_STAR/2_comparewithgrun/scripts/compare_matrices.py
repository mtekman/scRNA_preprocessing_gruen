#!/usr/bin/env python

import sys
from GeneralUtils import GeneralUtils as gu

if len(sys.argv) != 3:
    print('''
Compares two umi/barcode matrices of the same gene

%s <my_matrix.detailed> <grüns_matrix.detailed>

''' % sys.argv[0], file=sys.stderr)
    exit(-1)

#detailed = (len(sys.argv) == 4) and (sys.argv[3] == "--detailed")

mer = {} # barcode -> umi -> count or reads
umi_order = {}

def matrix2map(file, symbol):
    with open(file, 'r') as mat:
        string = ""
        maxline = 0
        while maxline < 100:
            line = mat.readline()
            if line.find('GENEID') != -1:
                break
            string += line
            maxline += 1

        headers = gu.header2barcode(string)

        for line in mat:
            tokens = line.split()
            gene, rbar = tokens[:2]

            umi_order[rbar] = True
            vals = tokens[2:]

            for i in range(len(headers)):
                bar = headers[i]
                val = vals[i]

                if bar not in mer:
                    mer[bar] = {}

                if rbar not in mer[bar]:
                    mer[bar][rbar] = {symbol : 0}

                mer[bar][rbar][symbol] = val


matrix2map(sys.argv[1],'M')
matrix2map(sys.argv[2],'G')

# Print matrix
barorder = list(mer.keys())
umiorder = list(umi_order.keys())

print("\nUmi(rows) vs Cells(cols):\n  overlap|reads(Mine),reads(Grün)\n")
print(' ' * 7, ' '.join(["%8s" % x for x in barorder]), sep="")

for umi in umiorder:
    print(umi, end="")
    for bar in barorder:
        out = None
        if umi in mer[bar]:
            res = mer[bar][umi]
            #if not(detailed):  # counts
            #    M = 0; G = 0
            #else:
            M = []; G = []; overlap = 0

            if 'M' in res:
                M = res['M'].split(';')
            if 'G' in res:
                G = res['G'].split(';')

            overlap = sum([x in M for x in G])

            #if umi == "GCTAGA" and bar == "ACTCTG":
            #    print(M, file=sys.stderr)
            #    print(G, file=sys.stderr)
            #    print(overlap, file=sys.stderr)
                
            out = "%2d:%2d,%2d" % (overlap, len(M),len(G))
        else:
            # Umi does not appear for that barcode for one file
            out = "--NONE--"
            
        print(" %s" % out, end="")
    print("")


#print(mer)
