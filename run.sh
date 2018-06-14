#! /bin/bash
#SBATCH -p bioinfo
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=20G
#SBATCH --export=ALL
#SBATCH --mail-user=sagar@ie-freiburg.mpg.de
#SBATCH --mail-type=ALL


f1="$1" #./WD_DP_120218_P1_1_R1.fastq"
f2="$2" # ./WD_DP_120218_P1_1_R2.fastq"
o="$3"  #WD_DP_120218_P1_1"
d="$3"  #WD_DP_120218_P1_1"

b="$4" #./celseq_barcodes.192.txt"
g="$5" #./Danio_rerio_Zv9_ens74_extended3_genes_ERCC92_GFP.fa"

[ "$5" = "" ] && echo "`basename $0` <R1.fastq> <R2.fastq> <output_dir> <barcodes.txt> <genome.fa>

The genome needs to be indexed, so run *bwa index -a bwtsw <genome.fa>* first

" && exit -1

t=4
i=0
rl=6
bl=6
te=0

./bin/do_mappings_strand.pl -test=$te -r=$g -f1=$f1 -f2=$f2 -out=$o -outdir=$d -t=$t -i=$i -cel=1 -fstr=1 -bar=$b -rb=1 -uniq=1 -rb_len=$rl >$o.map.log1 2>$o.map.log2
./bin/extract_counts_rb.pl -bl=$bl -in=$d/$o.cout.csv -outc=$d/$o.coutc.csv -outb=$d/$o.coutb.csv -outt=$d/$o.coutt.csv > $o.eco.log1 2> $o.eco.log2
