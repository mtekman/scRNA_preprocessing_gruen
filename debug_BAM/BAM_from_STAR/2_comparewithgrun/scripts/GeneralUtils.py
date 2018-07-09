#!/us/bin/env python

class GeneralUtils:
    umi_map = {}  # RBAR -> BAR -> num_reads
    barc_map = {} # BAR
    wanted_bars = {}

    @staticmethod
    def wantedBarcodes(barcode_file):
        with open(barcode_file, 'r') as wb:
            for line in wb:
                num, bar = line.split()
                GeneralUtils.wanted_bars[bar.splitlines()[0].strip()] = True

    @staticmethod
    def parseLoci(string):
        valid = True
        try:
            l_r = string.split(":")
            p1_p2 = l_r[1].split("-")

            chrom = l_r[0]
            loc1 = int(p1_p2[0])
            loc2 = int(p1_p2[1])

            if loc1 >= loc2:
                valid=False

        except (IndexError, ValueError) as e:
            valid = False

        if not(valid):
            print("Not a valid locus", string, file=sys.stderr)
            exit(-1)

        return(chrom,loc1,loc2)

    @staticmethod
    def updateMap(readname, umi=None, bar=None):
        gu = GeneralUtils
        purename = None

        if umi is None:
            tokens = readname.split("_")
            purename = tokens[0]
            bar = tokens[1]
            umi = tokens[2]
        else:
            purename = readname

        if bar not in gu.wanted_bars:
            return False

        if umi not in gu.umi_map:
            gu.umi_map[umi] = {}

        if bar not in gu.umi_map[umi]:
            gu.umi_map[umi][bar] = {}

        gu.barc_map[bar] = True
        gu.umi_map[umi][bar][purename] = True
        return True

    @staticmethod
    def printMap(gene_of_interest, outfilename):
        gu = GeneralUtils

        # Finally print out the map
        umi_order = list(gu.umi_map.keys())
        barc_order = list(gu.barc_map.keys())

        matrix_simple = open(outfilename +".matrix_simple", 'w')
        matrix_detail = open(outfilename +".matrix_detailed", 'w')

        header = "\n%23s%8s" % ("GENEID", "RBAR")
        barc_header =  gu.barcode2header(barc_order, 33)

        full_header = barc_header + header
        print(full_header, file=matrix_simple)
        print(full_header, file=matrix_detail)

        for umi in umi_order:
            out_nums = []
            out_read = []
            for bar in barc_order:
                try:
                    reads = gu.umi_map[umi][bar].keys()
                    readnames = ";".join(list(reads))
                    num = len(reads)  #
                except KeyError:
                    readnames = "."
                    num = 0

                out_nums.append(num)
                out_read.append(readnames)

            out_nums = '   '.join(["%2d" % o for o in out_nums])
            out_read = '   '.join(out_read)

            print("%23s%8s" % (gene_of_interest, umi), out_nums, sep=' ', file=matrix_simple)
            print("%23s%8s" % (gene_of_interest, umi), out_read, sep=' ', file=matrix_detail)

    @staticmethod
    def header2barcode(header):
        return("".join(
            ["".join(y) for y in list(zip(*[list(x) for x in header.splitlines()]))]
        ).split())

    @staticmethod
    def barcode2header(barcodes, buffer_size):
        buff = ' ' * buffer_size
        return('\n'.join(['%s%s' % (buff ,y) for y in [
            "    ".join(x) for x in list(zip(*[list(x) for x in barcodes]))]
        ]))
