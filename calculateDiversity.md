# Calculating Diversity at Synonymous and Nonsynonymous Positions, Including Invariable Sites
This markdown file organizes and comments on each step of the process to calculate diversity at synonymous and nonsynonymous positions, including invariable sites.

## Step 1: Filter Variants

```bash
vcftools --gzvcf /home/merrbii/Scratch/pprRelax/haplome/consensus_hetgt_soft_HI/geno_AB_snp_an12_ql20_star_hetgt_dp10.gvcf.gz --recode --recode-INFO-all --max-missing 0.8 --out hapAB.allPops.miss0.8 --remove-indels
```

Filters the VCF file to remove indels and keep variants with a maximum of 20% missing data.

## Step 2: Annotate Variants

```bash
java -jar ~/programs/snpEff/snpEff.jar -v -stats hapAB.allPops.miss0.8.ann.html Ppr hapAB.allPops.miss0.8.recode.vcf | bgzip > hapAB.allPops.miss0.8.ann.vcf.gz
```

Annotates the filtered VCF file using SnpEff and compresses the output.

## Step 3: Extract Invariable Sites

```bash
vcftools --gzvcf hapAB.allPops.miss0.8.ann.vcf.gz --recode --recode-INFO-all --max-maf 0 --out hapAB.allPops.miss0.8.ann.invar
```

Extracts invariable sites (sites with a minor allele frequency of 0).

## Step 4: Separate Synonymous and Nonsynonymous Variants

```bash
java -jar ~/programs/snpEff/SnpSift.jar filter "ANN[*].EFFECT has 'missense_variant'" hapAB.allPops.miss0.8.ann.vcf.gz > hapAB.allPops.miss0.8.ann.missense_any.vcf
java -jar ~/programs/snpEff/SnpSift.jar filter "ANN[*].EFFECT has 'synonymous_variant'" hapAB.allPops.miss0.8.ann.vcf.gz > hapAB.allPops.miss0.8.ann.synonymous_any.vcf
```

Filters the annotated VCF to separate synonymous and nonsynonymous variants.

## Step 5: Compress and Index VCF Files

```bash
ls *any.vcf | parallel --jobs 2 --eta bgzip {}
bgzip hapAB.allPops.miss0.8.ann.invar.recode.vcf
ls *vcf.gz | parallel --jobs 3 --eta tabix {}
```

Compresses and indexes the VCF files for further processing.

## Step 6: Concatenate VCF Files

```bash
bcftools concat --allow-overlaps hapAB.allPops.miss0.8.ann.invar.recode.vcf.gz hapAB.allPops.miss0.8.ann.synonymous_any.vcf.gz | bcftools view -e 'GT="."' -Oz -o hapAB.allPops.miss0.8.invar.synonymous_any.vcf.gz
bcftools concat --allow-overlaps hapAB.allPops.miss0.8.ann.invar.recode.vcf.gz hapAB.allPops.miss0.8.ann.missense_any.vcf.gz | bcftools view -e 'GT="."' -Oz -o hapAB.allPops.miss0.8.invar.missense_any.vcf.gz
```

Concatenates the invariable sites with the synonymous and nonsynonymous variants separately.

## Step 7: Calculate Nucleotide Diversity at Synonymous Sites

```bash
cut -f1 /RAID/Data/Mites/Genomes/Ppr/German_eiffel/hap0/Ppr.hap0.softmasked.fasta.fai > chr.txt

nice parallel -a chr.txt --colsep '\t' --jobs 10 python3 ~/programs/genomics_general/VCF_processing/parseVCF.py -i hapAB.allPops.miss0.8.invar.synonymous_any.vcf.gz --include {} --skipIndels -o geno/hapAB.allPops.miss0.8.invar.synonymous_any.{}.geno.gz

nice parallel -a chr.txt --colsep '\t' --jobs 10 python3 ~/programs/genomics_general/popgenWindows.py --windType sites -w 25000 -m 10000 --roundTo 6 -g geno/hapAB.allPops.miss0.8.invar.synonymous_any.{}.geno.gz -o csv/hapAB.allPops.miss0.8.invar.synonymous_any.{}.w25ksitesMin10k.csv.gz -f phased -T 4 -p CA -p CB -p GA -p GB -p IA -p IB -p JA -p JB -p RA -p RB --popsFile pops.txt --writeFailedWindows
```

Calculates nucleotide diversity at synonymous sites using a sliding window approach.
`pops.txt` is a two-column file with sampleID in the 1st column and population in the 2nd column. Please see Simon Martin's repository here: [https://github.com/simonhmartin/genomics_general](https://github.com/simonhmartin/genomics_general)

## Step 8: Calculate Nucleotide Diversity at Nonsynonymous Sites

```bash
nice parallel -a chr.txt --colsep '\t' --jobs 10 python3 ~/programs/genomics_general/VCF_processing/parseVCF.py -i hapAB.allPops.miss0.8.invar.missense_any.vcf.gz --include {} --skipIndels -o geno/hapAB.allPops.miss0.8.invar.missense_any.{}.geno.gz

nice parallel -a chr.txt --colsep '\t' --jobs 10 python3 ~/programs/genomics_general/popgenWindows.py --windType sites -w 25000 -m 10000 --roundTo 6 -g geno/hapAB.allPops.miss0.8.invar.missense_any.{}.geno.gz -o csv/hapAB.allPops.miss0.8.invar.missense_any.{}.w25ksitesMin10k.csv.gz -f phased -T 4 -p CA -p CB -p GA -p GB -p IA -p IB -p JA -p JB -p RA -p RB --popsFile pops.txt --writeFailedWindows
```

Calculates nucleotide diversity at nonsynonymous sites using a sliding window approach.

## Step 9: Calculate Overall Nucleotide Diversity

```bash
bcftools view -e 'GT="."' hapAB.allPops.miss0.8.ann.vcf.gz -Oz -o hapAB.allPops.miss0.8.ann.fixedGT.vcf.gz

nice parallel -a chr.txt --colsep '\t' --jobs 10 python3 ~/programs/genomics_general/VCF_processing/parseVCF.py -i hapAB.allPops.miss0.8.ann.fixedGT.vcf.gz --include {} --skipIndels -o geno/hapAB.allPops.miss0.8.ann.{}.geno.gz

nice parallel -a chr.txt --colsep '\t' --jobs 10 python3 ~/programs/genomics_general/popgenWindows.py --windType sites -w 25000 -m 10000 --roundTo 6 -g geno/hapAB.allPops.miss0.8.ann.{}.geno.gz -o csv/hapAB.allPops.miss0.8.ann.{}.w25ksitesMin10k.csv.gz -f phased -T 4 -p CA -p CB -p GA -p GB -p IA -p IB -p JA -p JB -p RA -p RB --popsFile pops.txt --writeFailedWindows
```

Calculates overall nucleotide diversity using a sliding window approach.
