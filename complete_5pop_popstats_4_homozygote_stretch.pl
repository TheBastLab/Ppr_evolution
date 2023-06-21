#!perl -w
use strict;

my $ind = $ARGV[0];
my $indx = (substr $ind,3,4)."_".(substr $ind,0,2);

my @b; my @minfo; my @cinfo; my $line;
my $mtemp; my $ctemp;
my $ix; my $chr; my $covnum;
my $nuc_mem; my $hom_leng;
my @out; my @result;

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

$nuc_mem = 1;
while ($line = <VCF>) {
 @b = split "\t", $line;
 if ($b[0] =~ /chr(\d+)/) {$chr = $1;}
 else {next;}

 if ((substr $b[$ix], 0, 1) ne ".") {
  if ((substr $b[$ix], 0, 1) != (substr $b[$ix], 2, 1)) {
   $hom_leng = $b[1] - $nuc_mem;
   if (($hom_leng > 0) & ($hom_leng <= 2000)) {
    $ctemp = substr $cinfo[$chr], ($nuc_mem-1), ($b[1]-$nuc_mem+1);
    $mtemp = substr $minfo[$chr], ($nuc_mem-1), ($b[1]-$nuc_mem+1);

    if ($ctemp =~ /^1+$/) {$out[0][$hom_leng] += 1;}
    elsif ($mtemp =~ /^1+$/) {
     $covnum = () = $ctemp =~ /1/g;
     if ($covnum > ($hom_leng/2)) {$out[1][$hom_leng] += 1;}
     else {$out[2][$hom_leng] += 1;}
    }
    else {$out[3][$hom_leng] += 1;}
    
   }
   
   
   $nuc_mem = $b[1];
  
 }}
 
}

push @result, "DbCHet\tFull_HighCov\tMean_HighCov\tOther_Nonrepeats\tRepeats\n";
foreach $hom_leng (1..2000) {
 push @result, $hom_leng;
 foreach (0..3) {
  if (defined $out[$_][$hom_leng]) {push @result, "\t$out[$_][$hom_leng]";}
  else {push @result, "\t0";}
 }
 push @result, "\n";
 
}

 open(FILEHANDLE, ">New_popstats/Homozygosity_stretch_stats_$ind.txt");
 print FILEHANDLE @result;

