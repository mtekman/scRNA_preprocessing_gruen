#!/usr/bin/perl -s -w

use lib './bin/';
use tools;

if (scalar @ARGV == 1)
{
    die "usage: -in=INPUTFILE.sam -s=1 if single end 0 is paired end -cmp=LIST_OF_READS.csv -s_flag= 1 or 0 ( if 1 then produce separate files for sense and antisense strand ) -u= 1 or 0 ( if 1 then only map reads that map to only one strand, optional) -gff=TRANSCRIPTOME.gff -uniq=1 (optional, only unique reads) -edit= maximum edit distance\n" if $ARGV[0] eq "help";
}

$gff = 0 if !$gff;
# gff initialization
$parent = "transcript";
$child  = "exon";
$check{$parent} = "ID";
$check{$child}  = "Parent";

$u = 0 if !$u;
$uniq = 0 if !$uniq;
$s_flag = 0 if !$s_flag;
$edit = 2 if !$edit;

%t = (); # hash for reads per reference sequence
%w = (); # hash for wig file
$n = 0; # number of reference sequences
$sh_pair = 0; # number of reads in pairs discarded due to length constraint
$sh_single = 0; # number of single reads in a pair discarded due to length constraint
$min_l = 15; # minimal length required for a read to be included
$u_r = 0; # number of unpaired reads
$p_r = 0; # number of paired reads
$u_m = 0; # number of reads mapped in a proper pair
$p_m = 0; # number of reads not mapped in a proper pair
$shared = 0; # number of reads shared with reference list in cmp if given
%seen = (); # hash for reads seen in reference list

@F = split(/\//,$in);
$in_name = $F[$#F];
$dir = join("/",@F[0..$#F - 1]);
$in_name =~ s/(\.)\w+$//;

$sout = $in; # output file for summary stats
if ($u){
    $sout =~ s/(\.)\w+$/\.u\.sout/;
}else{
    $sout =~ s/(\.)\w+$/\.sout/;
}
$rout = $in; # output file with list of mapped reads
if ($u){
    $rout =~ s/(\.)\w+$/\.u\.rout\.csv/;
}else{
    $rout =~ s/(\.)\w+$/\.rout\.csv/;
}
if ( $s_flag ){
    $cout_s = $in; # output file for read counts per reference sequence in input file
    $cout_a = $in; # output file for read counts per reference sequence in input file
    if ($u){
	$cout_s =~ s/(\.)\w+$/\.u\.sense\.cout\.csv/;
	$cout_a =~ s/(\.)\w+$/\.u\.antisense\.cout\.csv/;
    }else{
    	$cout_s =~ s/(\.)\w+$/\.sense\.cout\.csv/;
	$cout_a =~ s/(\.)\w+$/\.antisense\.cout\.csv/;
    }
}else{
    $cout   = $in;
    if ($u){
	$cout   =~ s/(\.)\w+$/\.u\.cout\.csv/;
    }else{
    	$cout   =~ s/(\.)\w+$/\.cout\.csv/;
    }
}
if ($cmp){
    open(CMP,'<',$cmp);
    while(<CMP>){
	chomp;
	$seen{$_} = 1;
    }
    close(CMP);
}

open(IN,'<',$in);
open(ROUT,'>',$rout);


$l_flag = 0;
$l_flag = 1 if $s;
while(<IN>){
    if (/^\@/){
      if (/^\@SQ/){
	chomp;
	@r = split(/\t/);
	$r[1] =~ s/SN://g;
	$t{"+"}{$r[1]} = 0;
	$t{"-"}{$r[1]} = 0;
	$n ++;
      }
    }else{
	if ($l_flag == 0){
	    @l = ($_);
	    $l_flag = 1;
	    next;
	}else{
	    if($s){
		@l = ($_);
	    }else{
		push(@l,$_);
		$l_flag = 0;
	    }
	    @d_flag = ();
	    for $i (0..$#l){
		($QNAME[$i],$FLAG[$i],$RNAME[$i],$POS[$i],$MAPQ[$i],$CIGAR[$i],
		 $MRNM[$i],$MPOS[$i],$ISIZE[$i],$SEQ[$i],$QUAL[$i])=split(/\t/,$l[$i]);
		@LINE = split(/\t/,$l[$i]);
		$NM[$i] = "NA";
		$XA[$i] = "NA";
		$d_flag[$i] = 1;
		foreach $el (@LINE){
		    ($dum,$dum,$NM[$i]) = split(/\:/,$el) if ($el =~ /^NM\:/);
		    ($dum,$dum,$XA[$i]) = split(/\:/,$el) if ($el =~ /^XA\:/);
		}
		##$NM[$i] = nb_mm($CIGAR[$i]) if $NM[$i] ne "NA";
		$d_flag[$i] = 0 if $XA[$i] eq "NA";
	    }
	    for $i (0..$#l){
		#next if $uniq && $d_flag[$i] == 1;
		if ( $uniq ){
		  @ALT  = split(/;/,$XA[$i]);
		  $u_flag = 0;
		  foreach $xa (@ALT){
		    last if $xa eq "NA" || $xa =~ /\s+/;
		    ($alt_rname, $pos, $CIGAR, $nm) = split(/,/,$xa);
		    $nm = nb_mm($CIGAR);
		    $u_flag = 1 if ($nm <= $NM[$i]);
		  }
		  next if $u_flag;
		}
		@flag   = split(//,reverse(dec2bin($FLAG[$i])));
		if ($flag[4]){$str = "-";}else{$str = "+";}
#		print STDERR join("",@flag)."\n";

#		if ($str eq "+"){
		$start  = $POS[$i] - 1;
		$length = length($SEQ[$i]);
#		}else{
#		    $length = length($SEQ[$i]);
#		    $start  = ( $POS[$i] - 1 ) - $length + 1;
#		}
#		print STDERR  join("",@flag)."\t".$str."\t".$RNAME[$i]."\t".$CIGAR[$i]."\t".$SEQ[$i]."\t".$start."\t".$length."\n";
		if (!$s){
		    if ($flag[0]){
			$p_r ++;
			if (length($SEQ[0]) < $min_l && length($SEQ[1]) < $min_l){
			    $sh_pair += 1;
			    next;
			}
			if (length($SEQ[$i]) < $min_l){
			    $sh_single += 1;
			    next;
			}
			if (length($SEQ[0]) < $min_l || length($SEQ[1]) < $min_l){
			    if ($flag[2] == 0 && $NM[$i] <= $edit){
				$u_m ++;
				update_t($RNAME[$i], $NM[$i], $XA[$i], \%t, $str, $u, $start, $length, \%w, $s_flag);
				$shared ++ if exists($seen{$QNAME[$i]});
				print ROUT $QNAME[$i]."\n";
				next;
			    }
			}
			if ($flag[2] == 0 && $NM[$i] <= $edit){
			    if ($flag[1]){
				$p_m ++;
			    }else{
				$u_m ++;
				
			    }
			    update_t($RNAME[$i], $NM[$i], $XA[$i], \%t, $str, $u, $start, $length, \%w, $s_flag, $gff);
			    $shared ++ if exists($seen{$QNAME[$i]});
			    print ROUT $QNAME[$i]."\n";
			}
		    }else{
			$u_r ++;
			next;
		    }
		}else{
		    $u_r ++;
		    if (length($SEQ[$i]) < $min_l){
			$sh_single += 1;
			next;
		    }
		    if ($flag[2] == 0 && $NM[$i] <= $edit){
			$u_m ++;
			update_t($RNAME[$i], $NM[$i], $XA[$i], \%t, $str, $u, $start, $length, \%w, $s_flag, $gff);
			$shared ++ if exists($seen{$QNAME[$i]});
			print ROUT $QNAME[$i]."\n";
		    }
		}
	    }
	}
    }
}
close(IN);

open(SOUT,'>',$sout);

print SOUT "number of reference sequences:\t".$n."\n";
$m = $p_r + $u_r;
print SOUT "number of reads:\t".$m."\n";
print SOUT "number of reads present in pairs:\t".$p_r ."\n";
print SOUT "number of reads not present in pairs:\t".$u_r."\n";
print SOUT "applied constraint on read length:\t".$min_l."\n";
print SOUT "number of reads in pairs discarded due to length constraint:\t".$sh_pair."\n";
print SOUT "number of single reads discarded due to length constraint:\t".$sh_single."\n";
$m = $p_m + $u_m;
print SOUT "number of mapped reads:\t".$m."\n";
print SOUT "number of reads mapped in pairs:\t".$p_m."\n";
print SOUT "number of reads mapped unpaired:\t".$u_m."\n";
if (!$s){
    $m = $p_m/$p_r;
    print SOUT "fraction of reads mapped in pairs:\t".$m."\n";
    $m = $u_m/$p_r;
    print SOUT "fraction of reads mapped unpaired:\t".$m."\n";
}else{
    $m = $u_m/$u_r;
    print SOUT "fraction of reads mapped unpaired:\t".$m."\n";
}
print SOUT "number of mapped reads present in reference sample:\t".$shared."\n";
$m = $shared/($p_m + $u_m);
print SOUT "fraction of mapped reads present in reference sample:\t".$m."\n";

$a_sum = 0;
$s_sum = 0;

if ($s_flag){
    open(COUTA,'>',$cout_a);
    open(COUTS,'>',$cout_s);
}else{
    open(COUT,'>',$cout);
}

foreach $k (keys %{$t{"+"}}){
    if ($s_flag){
	print COUTS $k."\t".$t{"+"}{$k}."\n";
	print COUTA $k."\t".$t{"-"}{$k}."\n";
    }else{
	print COUT $k."\t".($t{"-"}{$k} + $t{"+"}{$k})."\n";
    }
    $s_sum += $t{"+"}{$k};
    $a_sum += $t{"-"}{$k};
}
print SOUT "number of reads mapped to sense strand:\t".$s_sum."\n";
print SOUT "number of reads mapped to antisense strand:\t".$a_sum."\n";
if ( $s_sum + $a_sum > 0 ){
    $fr = $s_sum/($s_sum + $a_sum);
}else{
    $fr = 0;
}
print SOUT "fraction of reads mapped to sense strand:\t".$fr."\n";
if ( $s_sum + $a_sum > 0 ){
    $fr = $a_sum/($s_sum + $a_sum);
}else{
    $fr = 0;
}
print SOUT "fraction of reads mapped to antisense strand:\t".$fr."\n";

close(SOUT);
if ($s_flag){
    close(COUTA);
    close(COUTS);
}else{
    close(COUT);
}

if ($gff){
    open(IN,"<",$gff);
    while(<IN>){
	chomp;
	next if (/^\#\#gff/);
	last if (/^\#\#FASTA/);
	@F  = split(/\t/);
	next if $F[2] ne $child;
	$id = get_name_gff($_,\%check); 
	$exon{$id}{$F[3]} = $F[4];
	$chr{$id} = $F[0];
	$tr_str{$id} = $F[6];
    } 
    close(IN);
    for $id (keys %exon){
	$sum = 0;
	for $sc (sort {$a <=> $b} keys %{$exon{$id}}){
	    $shift{$id}{$sum} = $sc;
	    $sum += $exon{$id}{$sc} - $sc + 1; # results in greater or equal condition
	}
	$tr_length{$id} = $sum;
    }
    foreach $st (keys %w){
	foreach $id (sort keys %{$w{$st}}){
	    @ex_a = sort {$a <=> $b} keys %{$w{$st}{$id}};
	    if ($tr_str{$id} eq "-"){@ex_a = reverse(@ex_a);}
	    foreach $sct (@ex_a){
		$c_sh = 0;
		$sct_f = $sct;
		$sct_f = $tr_length{$id} - 1 - $sct if ($tr_str{$id} eq "-");
		foreach $sce (sort {$a <=> $b} keys %{$shift{$id}}){
#		    print $tr_str{$id}."\t".$id."\t".$sce."\t".$sct."\t".$sct_f."\t".$tr_length{$id}."\t".$shift{$id}{$sce}."\n";
		    if ($sct_f >= $sce){
			$c_sh = $shift{$id}{$sce} + $sct_f - $sce;
		    }else{
			last;
		    }
		}
#		print $tr_str{$id}."\t".$id."\t".$c_sh."\t".$sct."\t".$tr_length{$id}."\n";
#		if ( $c_sh == 0){
#		    print "**********\n";
#		    print join("\t",($st,$id,$chr{$id},$c_sh,$sct,$sct_f,$w{$st}{$id}{$sct}))."\n";
#		    foreach $sce (sort {$a <=> $b} keys %{$shift{$id}}){
#			print $id."\t".$sce."\t".$sct."\t".$sct_f."\t".$shift{$id}{$sce}."\n";
#		    }
#		    print "**********\n\n";
#		}
		if (!exists($w_f{$st}{$chr{$id}}{$c_sh})){
		    $w_f{$st}{$chr{$id}}{$c_sh}  = $w{$st}{$id}{$sct};
		}else{
		    $w_f{$st}{$chr{$id}}{$c_sh} += $w{$st}{$id}{$sct};
		}
#		print join("\t",($st,$id,$chr{$id},$c_sh,$w{$st}{$id}{$sct}))."\n";
	    }
	}
    }
    @F    = ("+", "-");
    @G    = ("sense", "antisense");
    if ($u) { $in_name = $in_name.".u";}
    @wout = ($in_name.".".$G[0].".wig",  $in_name.".".$G[1].".wig");
    foreach $i (0..$#F){
	$title = "track type=wiggle_0 name=\"".$in_name." (".$G[$i].")"."\" description=\"Read Profile ".$in_name." (".$G[$i].")"."\" visibility=full color=0,0,0 altColor=128,128,128 priority=N autoScale=on alwaysZero=off gridDefault=off maxHeightPixels=128:128:11 graphType=bar yLineMark=0.0 yLineOnOff=off windowingFunction=maximum smoothingWindow=off";
	if ( $dir ){
	    open(OUT,">",$dir."/".$wout[$i]);
	}else{
	    open(OUT,">",$wout[$i]);
	}
	print OUT $title."\n";
	foreach $c (sort keys %{$w_f{$F[$i]}}){
	    print OUT "variableStep chrom=".$c." span=1"."\n";
#	    print OUT "variableStep chrom=chr".$c." span=1"."\n";
	    foreach $st (sort {$a <=> $b} keys %{$w_f{$F[$i]}{$c}}){
		next if $st == 0;
		#print join("\t",($c,$st,$w_f{$F[$i]}{$c}{$st}))."\n";
		print OUT $st." ".$w_f{$F[$i]}{$c}{$st}."\n";
	    }
	}
	close(OUT);
    }
}

sub nb_mm {
  my $CIGAR = shift;
  return(0) if ! /S|I|D/;
  @x = split(//,$CIGAR);
  $k = $x[0];
  $mm = 0;
  for $i (1..$#x){
    if ( $x[$i] =~ /\d/ ){
      $k = $k.$x[$i];
    }elsif ( $x[$i] eq "M" ){
      $k = "";
    }elsif ( $x[$i] eq "S" ){
      $mm += $k;
      $k = "";
    }elsif ( $x[$i] eq "D" ){
      $mm += 1;
      $k = "";
    }elsif ( $x[$i] eq "I" ){
      $mm += 1;
      $k = "";
    }else{
      die "unknown symbol in CIGAR: ".$CIGAR."\n";
    }
  }
  return($mm);
}

sub update_t {
    my $RNAME = shift;
    my $NM    = shift;
    my $XA    = shift;
    my $t     = shift;
    my $st    = shift;
    my $u     = shift;
    my $start = shift;
    my $length= shift;
    my $w     = shift;
    my $flag  = shift;
    my $gff   = shift;
    my @ALT   = split(/;/,$XA);
    my %HITS   = ();
    my %COORDS = ();
    my $h_nb  = 1;
    fill_histo(\%{$HITS{$st}},$RNAME);

    foreach $xa (@ALT){
	last if $xa eq "NA" || $xa =~ /\s+/;
	my ($alt_rname, $pos, $CIGAR, $nm) = split(/,/,$xa);
	$nm = nb_mm($CIGAR);
	if ($pos < 0){$xa_st = "-";}else{$xa_st = "+";}
	if ($nm <= $NM){
	    fill_histo(\%{$HITS{$xa_st}},$alt_rname);
	    $h_nb ++; 
	}
    }
    if (!$u || !(exists($HITS{"+"}) && exists($HITS{"-"}))){
	foreach $s (keys %HITS){
	    foreach $rname (keys %{$HITS{$s}}){
		$$t{$s}{$rname} += $HITS{$s}{$rname}/$h_nb;				    
	    }
	}
	if ($gff){
	    for $i ($start .. $start + $length - 1){
		$st = "+" if (!$flag);
#		print "***".$RNAME."\t".$start."\t".$length."\n" if $i == 1364;
		if ( !exists($$w{$st}{$RNAME}{$i})){
		    $$w{$st}{$RNAME}{$i} = 1/$h_nb;
		}else{
		    $$w{$st}{$RNAME}{$i} += 1/$h_nb;
		}
	    }
	    foreach $xa (@ALT){
		last if $xa eq "NA" || $xa =~ /\s+/;
		my ($alt_rname, $pos, $CIGAR, $nm) = split(/,/,$xa);
		$nm = nb_mm($CIGAR);
		if ($pos < 0){$xa_st = "-";}else{$xa_st = "+";}
		$pos < 0 ? ($pos = -$pos) : ($pos = $pos);
		my $xa_start = $pos - 1;
		my $xa_stop  = $xa_start + $length - 1;
		if ($xa_st eq "-"){
#		    print STDERR "here ".$alt_rname." ".($pos - 1)."\n";
		    my $xa_stop  = $pos - 1;
		    my $xa_start = $xa_stop - $length + 1;
		}
		if ($nm <= $NM){
		    $xa_st = "+" if (!$flag);
#		    print STDERR $xa_start." ".$xa_stop."\n";
		    for $i ($xa_start .. $xa_stop){
#			print "***".$alt_rname.$xa_st."\t"."\t".$xa_start."\t".$xa_stop."\n" if $i == 1364;
			if (!exists($$w{$xa_st}{$alt_rname}{$i})){
			    $$w{$xa_st}{$alt_rname}{$i} = 1/$h_nb;
			}else{
			    $$w{$xa_st}{$alt_rname}{$i} += 1/$h_nb;
			}
		    }
		}
	    }
	}
    }
}

