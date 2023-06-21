#!perl -w
use strict;

my $ind1 = $ARGV[0];
my $ind2 = $ARGV[1];
my $indx1 = (substr $ind1,3,4)."_".(substr $ind1,0,2);
my $indx2 = (substr $ind2,3,4)."_".(substr $ind2,0,2);

my @b; my @cinfo1; my @cinfo2; my @cinfo;
my $line;
my $ix1; my $ix2; my $temp;
my $chr; my $blk; my $chr_blk; my $cst; my $mst;
my @dot; my @het; my $leng;
my @result;

open (FASTA, "Coverage/$ind1\_coverage_info.fa") or die; $chr = 0;
while ($line = <FASTA>) {
 if ($line =~ /^>/) {$chr += 1; next;}
 chomp $line;
 $cinfo1[$chr] .= $line;
}

open (FASTA, "Coverage/$ind2\_coverage_info.fa") or die; $chr = 0;
while ($line = <FASTA>) {
 if ($line =~ /^>/) {$chr += 1; next;}
 chomp $line;
 $cinfo2[$chr] .= $line;
}

foreach $chr (1..9) {
 foreach (1..(length ($cinfo1[$chr]))) {
  if ((substr $cinfo1[$chr], ($_-1), 1) + (substr $cinfo2[$chr], ($_-1), 1) == 2) {$cinfo[$chr] .= "1";}
  else {$cinfo[$chr] .= "0";}
 }

}

open (VCF, "New_popstats/Ppr_gatk_all25.snp.vcf");

foreach (1..45) {$line = <VCF>;}

$line = <VCF>; chomp $line; @b = split "\t", $line;
foreach (0..$#b) {
 if ($b[$_] eq $indx1) {$ix1 = $_;}
 if ($b[$_] eq $indx2) {$ix2 = $_;}
}

unless (defined $ix1) {die;}
unless (defined $ix2) {die;}

while ($line = <VCF>) {
 @b = split "\t", $line;
 if ($b[0] =~ /chr(\d+)/) {$chr = $1;}
 else {next;}

 $blk = int($b[1]/1000000);
 $chr_blk = $chr*100+$blk;
 
 if ((substr $cinfo[$chr], ($b[1]-1), 1) == 0) {next;}
 
 if (((substr $b[$ix1], 0, 1) eq ".") or ((substr $b[$ix2], 0, 1) eq ".")) {
  $dot[$chr_blk] += 1;
 }
 elsif ((substr $b[$ix1], 0, 1) != (substr $b[$ix1], 2, 1)) {
  if ((substr $b[$ix2], 0, 1) != (substr $b[$ix2], 2, 1)) {$het[0][$chr_blk] += 1;} #Het/Het
  else {$het[1][$chr_blk] += 1;} #Het/Hom
 }
 else {
  if ((substr $b[$ix2], 0, 1) != (substr $b[$ix2], 2, 1)) {$het[2][$chr_blk] += 1;} #Hom/Het
  elsif ((substr $b[$ix1], 0, 1) != (substr $b[$ix2], 0, 1)) {$het[3][$chr_blk] += 1;} #Hom/Hom with different alleles
 }

}

push @result, "Chr\tBlock\tLength_\tDots\tHet_Het\tHet_Hom\tHom_Het\tHom_dHom\n";

foreach $chr_blk (0..999) {
 unless (defined $het[0][$chr_blk]) {next;}
 
 $chr = int($chr_blk/100);
 $blk = $chr_blk % 100;

 if (1000000*($blk+1) > length($cinfo[$chr])) {$temp = substr $cinfo[$chr], (1000000*$blk);}
 else {$temp = substr $cinfo[$chr], (1000000*$blk), 1000000;}
 $leng = () = $temp =~ /1/g;

 push @result, "$chr\t$blk\t$leng";
 
 if (defined $dot[$chr_blk]) {push @result, "\t$dot[$chr_blk]";}
 else {push @result, "\t0";}

 if (defined $het[0][$chr_blk]) {push @result, "\t$het[0][$chr_blk]";}
 else {push @result, "\t0";}

 if (defined $het[1][$chr_blk]) {push @result, "\t$het[1][$chr_blk]";}
 else {push @result, "\t0";}

 if (defined $het[2][$chr_blk]) {push @result, "\t$het[2][$chr_blk]";}
 else {push @result, "\t0";}

 if (defined $het[3][$chr_blk]) {push @result, "\t$het[3][$chr_blk]";}
 else {push @result, "\t0";}

 push @result, "\n";
}

 open(FILEHANDLE, ">New_popstats/Shared_heterozygosity_$ind1\_$ind2.txt");
 print FILEHANDLE @result;
