#!perl -w
use strict;

my $ind = $ARGV[0];
my $thres = $ARGV[1]; #Coverage must be larger or equal to this number.
my $chr; my @b; my $nuc; my $line;
my $x1; my $x0; my $temp;
my @leng = (0,32608290,30492103,27587832,23397169,22442908,21557542,19673716,19511964,18765566);
my @genome;
my @result;
my @info; my @st;
my @minfo;

my $infile = "$ind\_coverage.bedgraph";
my $outfile = "$ind\_coverage_info.fa";

open (FASTA, "../Ppr_DE_Assembly_softmask.fa"); $chr = 0;
while ($line = <FASTA>) {
 if ($line =~ /^>/) {$chr += 1; next;}
 
 chomp $line;
 $line =~ s/[A-Z]/1/g;
 $line =~ s/[a-z]/0/g;
 $minfo[$chr] .= $line;
}
	print "mask info done\n";
open (FILE, $infile);

while ($line = <FILE>) {
 chomp $line; @b = split "\t", $line;
 if ($b[0] =~ /chr(\d+)/) {
  $chr = $1;
  if ($b[3] < $thres) {
   $temp = "0"x($b[2]-$b[1]);
   substr $minfo[$chr], $b[1], ($b[2]-$b[1]), $temp;
  }
 
 }
 
}

foreach (0..13) {$info[$_] = 0;}

foreach $chr (1..9) {

 @result = ();
 push @result, ">chr$chr\n";
 $temp = 0;
 while (($temp+60) < length($minfo[$chr])) {
  push @result, (substr $minfo[$chr], $temp, 60);
  push @result, "\n";
  $temp += 60;
 }
 push @result, (substr $minfo[$chr], $temp);
 push @result, "\n";
 open(FILEHANDLE, ">>$outfile");
 print FILEHANDLE @result;

}



