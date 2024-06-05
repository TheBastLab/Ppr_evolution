# Calculate evolutionary rate using haplo- population-sepecific variation and Hermannia gibba as outgroup sepecies:


## Step 1: Filtering Variants

```bash
vcftools --gzvcf /home/merrbii/Scratch/hyphy/consensus_hetgt_soft_HI/geno_AB_snp_an12_ql20_star_hetgt_dp10.gvcf.gz --recode --recode-INFO-all --keep samples.DEJA.txt --max-missing 0.9 --out samples.DEJA.miss0.9
```
The vcftools command filters variants to include only those present in samples from Japan and Germany with a maximum missing data threshold of 10%. The file samples.DEJA.txt contains samples from Japan and Germany.

## Step 2: Identify Fixed Alleles

```bash
for p in germany japan ; do 
    for h in A B; do 
        echo "bcftools view -S pops/${p}_${h}.txt samples.DEJA.miss0.9.recode.vcf | bcftools view -q 0.55 -Oz -o deja/${p}_${h}.DEJA.miss0.9.q0.55.vcf.gz";
    done;
done | parallel --jobs 4 --eta
```
This step identifies almost fixed alternate alleles in each haplotype/population. If not fixed for the alternate (ALT) allele, it must be fixed for the reference (REF) allele.

## Step 3: Basic Statistics

```bash
ls *gz | parallel --jobs 10 'tabix {}'

for i in germany japan ; do 
    echo "bcftools isec ${i}_A.DEJA.miss0.9.q0.55.vcf.gz ${i}_B.DEJA.miss0.9.q0.55.vcf.gz -p ${i}";
done | parallel --jobs 5 --eta
```

Indexes the VCF files and performs intersection analysis to identify shared and unique variants.

## Step 4: Obtain Consensus Sequences

```bash
for i in germany_A germany_B japan_A japan_B; do 
    echo "bcftools consensus -o hyphy/ref/${i}.fa -f /RAID/Data/Mites/Genomes/Ppr/German_eiffel/hap0/Ppr.hap0.softmasked.fasta ${i}.DEJA.miss0.9.q0.55.vcf.gz";
done | parallel --jobs 6 --eta

Obtains consensus sequences for each population.

## Step 5: Extract CDS

```bash
ls hga.fasta | nice parallel --jobs 10 --eta agat_sp_extract_sequences.pl -g ../gtf/{.}.gff3 -f {} -t cds -o ../cds/{.}.cds.fa

ls *B.fa *A.fa | nice parallel --jobs 10 --eta --dry-run agat_sp_extract_sequences.pl -g ../gtf/ppr.gff3 -f {} -t cds -o ../cds/{.}.cds.fa
```
Extracts coding sequences (CDS) from the reference genome.

## Step 6: Extract Protein Sequences

```bash
ls hga.fasta | nice parallel --jobs 10 --eta agat_sp_extract_sequences.pl -g ../gtf/{.}.gff3 -f {} -t cds -p -o ../pep/{.}.pep.fa

ls *B.fa *A.fa | nice parallel --jobs 10 --eta agat_sp_extract_sequences.pl -g ../gtf/ppr.gff3 -f {} -t cds -p -o ../pep/{.}.pep.fa
```
Extracts protein sequences from the CDS.

## Step 7: Get Longest Peptide

```bash
ls *fa | sed 's/.pep.reh.fa//g' > ../speciesList.txt

for i in $(cat ../speciesList.txt); do 
    echo "seqkit fx2tab -l ${i}.pep.reh.fa | sed 's/_gene:/\t/g' | sort -k2,2 -k4,4nr | sort -k2,2 -u -s | awk '{print \">\"\"\$1\";gene=\"\"\$2\";length=\"\"\$4\"\n\"\"\$3}' | sed 's/_indv:/;indv=/g' > longest.${i}.pep.reh.fa";
done > get.longest.iso.cmd

cat get.longest.iso.cmd | parallel --jobs 10 --eta
```
Gets the longest peptide sequences.

## Step 8: Get Longest CDS

```bash
for i in $(cat ../speciesList.txt); do 
    grep "^>" longest.$i.pep.reh.fa | cut -f1 -d";" | sed 's/>//g' > ../cds/$i.longest_cds.txt;
done

for i in *.fa; do 
    awk '/^>/ {$0=$1} 1' $i > tmp.fa && mv tmp.fa $i;
done

for i in $(cat ../speciesList.txt); do 
    seqkit grep -n -f $i.longest_cds.txt $i.cds.fa > longest.$i.cds.fa;
done
```
Gets the longest coding sequences.

## Step 9: Orthofinder

```bash
mkdir orthofinder && mkdir orthofinder/inputs

mv pep/longest.* orthofinder/inputs/

orthofinder -f inputs/ -t 40 > orthofinder.DEJA.log
```
Runs OrthoFinder to identify orthologous groups.

## Step 10: Prepare for Alignment

```bash
cat inputs/OrthoFinder/Results_May26/Orthogroups/Orthogroups_SingleCopyOrthologues.txt | while read group ; do 
    echo "cp inputs/OrthoFinder/Results_May26/Orthogroup_Sequences/${group}.fa alignments/${group}.faa";
done | parallel --jobs 40 --eta

cat inputs/OrthoFinder/Results_May26/Orthogroups/Orthogroups_SingleCopyOrthologues.txt | while read group ; do 
    echo "cat ./alignments/${group}.faa | grep '>' | sed 's/>//g' | seqkit grep -n -f - <(cat /home/merrbii/Scratch/hyphy/allPops/0.7Miss/sharedSNPs/vcfs/deja/hyphy/cds/*.fa) >> alignments/${group}.fna";
done > get.SCO.cds.cmd

cat get.SCO.cds.cmd | parallel --jobs 40 --eta
```
Prepares data for multiple sequence alignment.

## Step 11: Align

```bash
cat ../inputs/OrthoFinder/Results_May26/Orthogroups/Orthogroups_SingleCopyOrthologues.txt | parallel --jobs 70 --eta 'clustalo -i {1}.faa -o {1}.aln.faa'
```
Performs multiple sequence alignment.

## Step 12: Convert Alignment to Codon Alignment

```bash
cat ../inputs/OrthoFinder/Results_May26/Orthogroups/Orthogroups_SingleCopyOrthologues.txt | parallel --jobs 70 --eta pal2nal.pl {1}.aln.faa {1}.fna -output fasta -nogap ">" {1}.pal2nal.fasta "2>" {1}.pal2nal.fasta.log
```
Converts protein alignments to codon alignments.

## Step 13: Quality Control

```bash
grep "ERROR" *.log | sed 's/.pal2nal.*//g' > failed.txt

cat ../inputs/OrthoFinder/Results_May26/Orthogroups/Orthogroups_SingleCopyOrthologues.txt | grep -v -f failed.txt > passedAlignments.txt
```

Identifies errors in the alignment process and filters out failed alignments.

## Step 14: Phylogenetic Tree Construction

```bash
ls *.fasta | nice parallel --eta --jobs 40 'raxmlHPC -m GTRGAMMA -p 12345 -# 20 -s {} -n {.}.tree1 -T 1'

ls *.fasta | nice parallel --eta --jobs 40 'raxmlHPC -m GTRGAMMA -p 12345 -b 12345 -# 100 -s {} -n {.}.tree2 -T 1'

ls *.fasta | nice parallel --eta --jobs 40 'raxmlHPC -m GTRCAT -p 12345 -f b -t RAxML_bestTree.{.}.tree1 -z RAxML_bootstrap.{.}.tree2 -n {.}.tree3 -T 1'
```
Constructs phylogenetic trees using RAxML with different models and bootstrapping.

## Step 15: Clean up
```bash
for i in {0..19}; do 
    rm *.RUN.$i;
done

mkdir sequence && mv *.fasta* sequence

mkdir bestTrees && mkdir bestTrees/tmp

mv *.tree1 bestTrees/tmp/

mv *.tree2 bestTrees/tmp/

mv RAxML_info*.tree3 bestTrees/tmp/

mv RAxML_bipartitions.OG0* bestTrees/tmp/

mv RAxML_bipartitionsBranchLabels.OG0* bestTrees/
```

## Step 16: Fit Evolutionary Models

```bash
# Run fitmg94
mkdir FitMG94

cat passedAlignments.txt | while read group; do 
    echo "hyphy ~/programs/hyphy-analyses/FitMG94/FitMG94.bf --alignment sequence/${group}.pal2nal.fasta --tree bestTrees/RAxML_bipartitionsBranchLabels.${group}.pal2nal.tree3 --output FitMG94/fitMG94-hyphy-${group}.json --type local >> FitMG94/${group}.fitMG94.log 2>&1";
done | nice parallel --jobs 50 --eta
```
Fits the MG94 evolutionary model to each alignment using HyPhy.

## Step 17: Check for Errors!!

```bash 
grep "error" *.log | cut -f1 -d "." | sort | uniq
```
Checks for errors in the fitting process.
