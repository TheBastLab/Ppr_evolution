#!perl -w
use strict;

my $ind = $ARGV[0];
my $indx = (substr $ind,3,4)."_".(substr $ind,0,2);

my @b; my @minfo; my @cinfo; my $line;
my $ix; my $temp;
my $chr; my $blk; my $chr_blk; my $cst; my $mst;
my @dot; my @het; my @leng;
my @result;

open (FASTA, "Ppr_DE_Assembly_softmask.fa"); $chr = 0;
while ($line = <FASTA>) {
 if ($line =~ /^>/) {$chr += 1; next;}
 
 chomp $line;
 $line =~ s/[A-Z]/1/g;
 $line =~ s/[a-z]/0/g;
 $minfo[$chr] .= $line;
}


open (FASTA, "Coverage/$ind\_coverage_info.fa") or die; $chr = 0;
while ($line = <FASTA>) {
 if ($line =~ /^>/) {$chr += 1; next;}
 chomp $line;
 $cinfo[$chr] .= $line;
}


open (VCF, "New_popstats/Ppr_gatk_all25.snp.vcf");

foreach (1..45) {$line = <VCF>;}

$line = <VCF>; chomp $line; @b = split "\t", $line;
foreach (0..$#b) {if ($b[$_] eq $indx) {$ix = $_;}}

unless (defined $ix) {die;}

while ($line = <VCF>) {
 @b = split "\t", $line;
 if ($b[0] =~ /chr(\d+)/) {$chr = $1;}
 else {next;}

 $blk = int($b[1]/1000000);
 $chr_blk = $chr*100+$blk;
 
 $cst = substr $cinfo[$chr], ($b[1]-1), 1;
 $mst = substr $minfo[$chr], ($b[1]-1), 1;
 
 if ((substr $b[$ix], 0, 1) eq ".") {
  if ($cst == 1) {$dot[0][$chr_blk] += 1;}
  if ($mst == 1) {$dot[1][$chr_blk] += 1;}
  $dot[2][$chr_blk] += 1;
 }
 elsif ((substr $b[$ix], 0, 1) != (substr $b[$ix], 2, 1)) {
  if ($cst == 1) {$het[0][$chr_blk] += 1;}
  if ($mst == 1) {$het[1][$chr_blk] += 1;}
  $het[2][$chr_blk] += 1;
 }

}


push @result, "Chr\tBlock\tLength_HighCov\tHet_Highcov\tDots_Highcov\tLength_Mask\tHet_Mask\tDots_Mask\tLength_All\tHet_All\tDots_All\n";

foreach $chr_blk (0..999) {
 unless (defined $het[0][$chr_blk]) {next;}
 
 $chr = int($chr_blk/100);
 $blk = $chr_blk % 100;
 @leng = (0,0,0,0,0,0);
 
 if (1000000*($blk+1) > length($cinfo[$chr])) {$temp = substr $cinfo[$chr], (1000000*$blk);}
 else {$temp = substr $cinfo[$chr], (1000000*$blk), 1000000;}
 $leng[0] = () = $temp =~ /1/g;
 
 if (1000000*($blk+1) > length($cinfo[$chr])) {$temp = substr $minfo[$chr], (1000000*$blk);}
 else {$temp = substr $minfo[$chr], (1000000*$blk), 1000000;}
 $leng[1] = () = $temp =~ /1/g;
 
 if (1000000*($blk+1) > length($cinfo[$chr])) {$leng[2] = (length($cinfo[$chr]))- (1000000*$blk);}
 else {$leng[2] = 1000000;}
 
 

 push @result, "$chr\t$blk";
 foreach (0..2) {
  unless (defined ($het[$_][$chr_blk])) {$het[$_][$chr_blk] = 0;}
  unless (defined ($dot[$_][$chr_blk])) {$dot[$_][$chr_blk] = 0;}
  push @result, "\t$leng[$_]\t$het[$_][$chr_blk]\t$dot[$_][$chr_blk]";
 }
 push @result, "\n";

}

 open(FILEHANDLE, ">New_popstats/Heterozygosity_by_blocks_$ind.txt");
 print FILEHANDLE @result;

