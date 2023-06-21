#!perl -w
use strict;

my $ind = $ARGV[0];
my $chr; my @b; my $nuc; my $line;
my $x1; my $x0; my $temp;
my @cov_1; my @cov_0;
my @result;
my @info; my @st;
my @minfo;
my $cmax = 0;

my $infile = "$ind\_coverage.bedgraph";

open (FASTA, "../Ppr_DE_Assembly_softmask.fa"); $chr = 0;
while ($line = <FASTA>) {
 if ($line =~ /^>/) {$chr += 1; next;}
 
 chomp $line;
 $line =~ s/[A-Z]/1/g;
 $line =~ s/[a-z]/0/g;
 $minfo[$chr] .= $line;
}

open (FILE, $infile);

while ($line = <FILE>) {
 chomp $line; @b = split "\t", $line;
 if ($b[0] =~ /chr(\d+)/) {
  $chr = $1;
  if ($b[3] > $cmax) {$cmax = $b[3];}
  
  foreach $nuc ($b[1]..($b[2]-1)) {
   if ((substr $minfo[$chr], $nuc, 1) == 1) {$cov_1[$b[3]] += 1;}
   else {$cov_0[$b[3]] += 1;}
  }
 
 }
 
}

push @result, "Coverage\tNon_repeat\tRepeat\n";
foreach (0..$cmax) {
 unless ((defined $cov_1[$_]) or (defined $cov_0[$_])) {next;}

 push @result, "$_\t";
 if (defined $cov_1[$_]) {push @result, "$cov_1[$_]\t";}
 else {push @result, "0\t";}
 if (defined $cov_0[$_]) {push @result, "$cov_0[$_]\n";}
 else {push @result, "0\n";}

}


 open(FILEHANDLE, ">$ind\_coverage_distribution.txt");
 print FILEHANDLE @result;



