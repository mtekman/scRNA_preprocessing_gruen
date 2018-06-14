#!/usr/bin/perl -w -s


if (scalar @ARGV == 1 && $ARGV[0] eq "help"){
    die "usage: -in=INPUT.qseq_txt -clean=1 for 3' end trimming of base qualities B\n";
}
$clean=0 if !$clean;
$l = 0;
open(IN, $in);
while (<IN>){
    next if (substr($_, 0, 1) eq "#");
    chomp();
    @f = split(/\t/);
    warn "Problem at line ".$l.":".$_ if (scalar @f < 11);
    $f[8] =~ tr/\./N/;
    $l_2 = $f[8];
    $l_4 = $f[9];
    if ($clean){
	if ($l_4 =~ /B+$/){
	    $length = length($l_4) - length($&);
	    $length = 1 if $length == 0; 
	    $l_2 = substr($l_2,0,$length);
	    $l_4 = substr($l_4,0,$length);
	}
    }
    print "@" .$f[0].".".$f[1].".".$f[2].".".$f[3].".".$f[4].".".$f[5]."\n";
    print $l_2."\n";
    print "+\n";  
    print $l_4."\n";
}
close(IN);


