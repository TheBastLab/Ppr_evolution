### data needed for variant calling
##### reference genome 
Ppr.hap0.softmasked.fasta
##### bam files
Ppr*sorted.removed_duplicates.bam #generated with mapping.sh

## Prepare the data
#### build dict for genome
```
gatk-4.1.9.0/gatk CreateSequenceDictionary \
-R Ppr.hap0.softmasked.fasta
```
```
samtools faidx Ppr.hap0.softmasked.fasta
```
#### generate bam files 
```
mapping.sh:
```
##### mapping with bwa
```
bwa mem -t 40 -R "@RG\tID:${i}\tSM:${i}\tLB:${i}\tPL:Illumina" \
Ppr.hap0.softmasked.fasta \
R2_${i}_raw_val_2.fq.gz \
R1_${i}_raw_val_1.fq.gz > $i.sam
```
##### compress to bam and sort
```
samtools view -@ 40 -bS $i.sam > $i.bam
samtools sort -@ 40 $i.bam -o $i.sorted.bam
rm -f $i.sam
```
##### mark and remove duplication
```
java -jar picard.jar MarkDuplicates \
I=$i.sorted.bam \
O=$i.sorted.removed_duplicates.bam \
M=$i.removed_dup_metrics.txt \
REMOVE_DUPLICATES=true \
REMOVE_SEQUENCING_DUPLICATES=true
```
```
samtools index -@ 40 $i.sorted.removed_duplicates.bam
rm -f $i.sorted.bam
rm -f $i.bam
```
##### compare duplication removed bam with regular bam (not part of regular workflow)
```
comparebams.sh:
java -jar picard.jar CompareSAMs \
$i.sorted.bam \
$i.sorted.marked_duplicates.bam \
O=$i.comparison.tsv
```
##### get sample names
```
gatk-4.1.9.0/gatk GetSampleName \
     -I $i.sorted.bam \
     -O sample_name_$i.txt
```

# Variant calling

[GATK](https://gatk.broadinstitute.org/hc/en-us) version 4.1.9.0

In short: mapping population data to reference - HaplotypeCaller - CombineGVCF - GenotypeGVCFs - SelectVariants
##### Comment: Note down how many total variants were processed after each step.
##### Comment: run HaploypeCaller for each individual simultaneasly (to save time)
```sh
GATK HaplotypeCaller -R Ppr.hap0.softmasked.fasta \
        --emit-ref-confidence GVCF \
	--native-pair-hmm-threads 8 \
        -I $i.sorted.marked_duplicates.bam \
        -O Ppr_gatk_haplo$i
```
##### change header ID in vcf for vcf25 (see rename_vcf.sh)
```
java -jar /NVME/Software/picard.jar RenameSampleInVcf \
I=Ppr.RU.gatk_haploT507 \
O=Ppr.RU.gatk_haploT507.vcf \
OLD_SAMPLE_NAME=T507 \
NEW_SAMPLE_NAME=T507_RU
```
##### index
```
/NVME/Software/popgen/gatk-4.1.9.0/gatk IndexFeatureFile \
     -I *.vcf
```
##### merge 
```sh
GATK CombineGVCFs -R Ppr.hap0.softmasked.fasta \
	-V Ppr_gatk_haploT501 \
	-V Ppr_gatk_haploT502 \
	-V Ppr_gatk_haploT505 \
	-V Ppr_gatk_haploT506 \
	-V Ppr_gatk_haploT507 \
	-O Ppr_merged.g.vcf
```

```sh
GATK GenotypeGVCFs -R Ppr.hap0.softmasked.fasta \
        -V Ppr_merged.g.vcf \
        -O Ppr_gatk_all.vcf
```

```sh
GATK ValidateVariants --gvcf -R Ppr.hap0.softmasked.fasta \
	-V Ppr_gatk_all.vcf
```

```sh
GATK SelectVariants \
	-select-type SNP \
	-V Ppr_gatk_all.vcf.gz \
	-O Ppr_gatk_all.snp.vcf.gz
```

```sh
GATK VariantFiltration -V Ppr_gatk_all.snp.vcf.gz \
	--filter-expression "QD <2.0 || MQ <40.0 || FS >60.0 || SOR >3.0 || ReadPosRankSum < -8.0 || MQRankSum < -12.5" \
	--filter-name "PASS" \
	-O Ppr_gatk.SNP.filtered.vcf.gz
```

```sh
GATK VariantsToTable \
	-V Ppr_gatk.SNP.filtered.vcf.gz \
	-F CHROM -F POS -F REF -F ALT -F MULTI-ALLELIC -F TYPE -GF AD \
	-O Ppr_gatk.SNP.filtered.table
```

##### Counting the filtered SNPs
```
vcftools --gzvcf Ppr_gatk.SNP.filtered.vcf.gz --out Ppr_gatk.SNP.filtered.removed.vcf.gz --remove-filtered-all
```
