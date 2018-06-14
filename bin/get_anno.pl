#!/usr/bin/perl -s -w

use lib './bin/';
#use lib '/Users/d.grun/data/bin/';
use tools;

if (scalar @ARGV == 1)
{
    die "usage: -sam=INPUT.sam -st=1 (1: include strand information, 0: no strand information (default)" if $ARGV[0] eq "help";
}
$st=0 if !$st;

$l_flag = 0;
open(IN,"<",$sam);
while(<IN>){
  if (!/^\@/){
    if ($l_flag == 0){
      @l = ($_);
      $l_flag = 1;
      next;
    }else{
	$count ++;
      push(@l,$_);
      $l_flag = 0;
      @d_flag = ();
      for $i (0..$#l){
	($QNAME[$i],$FLAG[$i],$RNAME[$i],$POS[$i],$MAPQ[$i],$CIGAR[$i],
	 $MRNM[$i],$MPOS[$i],$ISIZE[$i],$SEQ[$i],$QUAL[$i])=split(/\t/,$l[$i]);
	@LINE = split(/\t/,$l[$i]);
	$NM[$i] = "NA";
	$XA[$i] = "NA";
	$d_flag[$i] = 1;
	foreach $el (@LINE){
	  ($dum,$dum,$NM[$i]) = split(/\:/,$el) if ($el =~ /^NM/);
	  ($dum,$dum,$XA[$i]) = split(/\:/,$el) if ($el =~ /^XA/);
	}
	$d_flag[$i] = 0 if $XA[$i] eq "NA";
      }
      for $i (0..$#l){
	@flag   = split(//,reverse(dec2bin($FLAG[$i])));
	if ($flag[4]){$str = "-";}else{$str = "+";}
	$start  = $POS[$i] - 1;
	$length = length($SEQ[$i]);
	update_cov($RNAME[$i], $NM[$i], $XA[$i], $str, $start, \%cov);
      }
    }
  }
}
close(IN);
      
sub update_cov {
  my $RNAME = shift;
  my $NM    = shift;
  my $XA    = shift;
  my $str   = shift;
  my $start = shift;
  my $cov   = shift;
  
  my @ALT   = split(/;/,$XA);
  $$cov{$RNAME}{$str}{$start} = 1;
  foreach $xa (@ALT){
    last if $xa eq "NA" || $xa =~ /\s+/;
    my ($alt_rname, $pos, $CIGAR, $nm) = split(/,/,$xa);
    if ($pos < 0){$xa_st = "-"; $pos = -$pos;}else{$xa_st = "+";}
    if ($nm <= $NM){
      $$cov{$alt_rname}{$xa_st}{$pos} = 1;
    }
  }
}


@anno = (
    "Mus_musculus.GRCm38.71.clean.snRNA.exon.gtf",
    "Mus_musculus.GRCm38.71.clean.snoRNA.exon.gtf",
    "mm10_tRNA.gtf",
    "Mus_musculus.GRCm38.71.clean.rRNA.exon.gtf",
    "mmu_prec_miRNA.gff3",
    "mmu_mat_miRNA.gff3",
    "Mus_musculus.GRCm38.71.clean.pseudogene.intron.gtf",
    "Mus_musculus.GRCm38.71.clean.pseudogene.exon.gtf",
    "Mus_musculus.GRCm38.71.clean.non_coding.intron.gtf",
    "Mus_musculus.GRCm38.71.clean.non_coding.exon.gtf",
    "Mus_musculus.GRCm38.71.clean.lincRNA.intron.gtf",
    "Mus_musculus.GRCm38.71.clean.lincRNA.exon.gtf",
    "Mus_musculus.GRCm38.71.clean.protein_coding.intron.gtf",
	 "Mus_musculus.GRCm38.71.clean.protein_coding.exon.gtf",
	 "mm10_RefSeq_genes_introns_clean.gtf",
	 "mm10_RefSeq_genes_clean_ERCC92.gtf"
    );


for $i (@anno) {
  open(IN,"<".$i);
  %seen =();
  while(<IN>){
    chomp;
    @F = split(/\s+/,$_);
    $name = get_entry_gtf($_,"gene_id") if ($i =~ /gtf/);
    $name = get_entry_gff($_,"Name"   ) if ($i =~ /gff/);

    for $j ($F[3] .. $F[4]){
      if ( $st ){
	next if !exists($cov{$F[0]}{$F[6]}{$j});
	%{$h{$F[0]}{$F[6]}{$j}} = () if !exists($seen{$F[0]}{$F[6]}{$j});
	$h{$F[0]}{$F[6]}{$j}{$name} = 1;
	$k{$F[0]}{$F[6]}{$j} = join(".",@F[1..2]) if $i =~ /gtf/;
	$k{$F[0]}{$F[6]}{$j} = $F[2] if $i =~ /gff/;
	$seen{$F[0]}{$F[6]}{$j} = 1;
      }else{
	next if !exists($cov{$F[0]}{"+"}{$j}) && !exists($cov{$F[0]}{"-"}{$j});
	for $str (("+","-")){
	  %{$h{$F[0]}{$str}{$j}} = () if !exists($seen{$F[0]}{$str}{$j});
	  $h{$F[0]}{$str}{$j}{$name} = 1;
	  $k{$F[0]}{$str}{$j} = join(".",@F[1..2]) if $i =~ /gtf/;
	  $k{$F[0]}{$str}{$j} = $F[2] if $i =~ /gff/;
	  $seen{$F[0]}{$str}{$j} = 1;
	}
      }
    }
  }
  close(IN);
}

foreach $k (sort keys %h){
    foreach $k1 (sort keys %{$h{$k}}){
	foreach $k2 (sort {$a <=> $b} keys %{$h{$k}{$k1}}){
	    print join("\t",($k,$k1,$k2,$k{$k}{$k1}{$k2},join(",",sort keys %{$h{$k}{$k1}{$k2}})))."\n";
	}
    }
}
