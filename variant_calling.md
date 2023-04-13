## software
gatk=/NVME/Software/popgen/gatk-4.1.9.0/gatk

## data
ref=/home/hoeztopr/Data/hoeztopr/Ppr/genomes/DE/Ppr.hap0.softmasked.fasta
bam=Ppr*sorted.removed_duplicates.bam #generated with mapping.sh


## Preparation
#### build dict for genome
```
/NVME/Software/popgen/gatk-4.1.9.0/gatk CreateSequenceDictionary \
-R /home/hoeztopr/Data/hoeztopr/Ppr/genomes/DE/Ppr.hap0.softmasked.fasta -O Ppr.dict
```
```
samtools faidx Ppr.italy.hap0.softmasked.fasta
```
#### generate bam files 
```
/home/hoeztopr/Data/hoeztopr/Scripts/mapping.sh:
```
#mapping (Italian Population as example)
```
bwa mem -t 40 -R "@RG\tID:${i}\tSM:${i}\tLB:${i}\tPL:Illumina" \
/RAID/Data/Mites/Genomes/Ppr/Italy/Ppr.italy.hap0.softmasked.fasta \
/RAID/Data/Mites/Reads/TELLSeq/Corrected_Reads/Ppr/Italy/German_ref_cleanReads/03Italy_R2_${i}_raw_val_2.fq.gz \
/RAID/Data/Mites/Reads/TELLSeq/Corrected_Reads/Ppr/Italy/German_ref_cleanReads/03Italy_R1_${i}_raw_val_1.fq.gz > $i.IT.sam
```
#compress to bam and sort
```
samtools view -@ 40 -bS $i.IT.sam > $i.IT.bam
samtools sort -@ 40 $i.IT.bam -o $i.IT.sorted.bam
rm -f $i.IT.sam
```
#mark and remove duplication
```
java -jar /NVME/Software/picard.jar MarkDuplicates \
I=$i.IT.sorted.bam \
O=$i.IT.sorted.removed_duplicates.bam \
M=$i.IT.removed_dup_metrics.txt \
REMOVE_DUPLICATES=true \
REMOVE_SEQUENCING_DUPLICATES=true
```
```
samtools index -@ 40 $i.IT.sorted.removed_duplicates.bam
rm -f $i.IT.sorted.bam
rm -f $i.IT.bam
```
##### compare duplication removed bam with regular bam 
(not part of regular workflow)
	```
	/home/hoeztopr/Data/hoeztopr/Scripts/comparebams.sh:
        java -jar /NVME/Software/picard.jar CompareSAMs \
	$i.sorted.bam \
	$i.sorted.marked_duplicates.bam \
	O=$i.comparison.tsv
	```
	
### get sample names (not necessary for workflow)
```
/NVME/Software/popgen/gatk-4.1.9.0/gatk GetSampleName \
     -I /RAID/Data/gaoshan/gaoshan/hifiasm_tell-seq/tell_sort/tell_sort_out/Ppr_T502/Ppr_T502_temp/Ppr_T502.sorted.bam \
     -O sample_nameT502.txt
```
### check if bam is good to go for GATK 
you can also skip this
```
gatk ValidateSamFile --INPUT=INDIVIDUAL_readgroup_fixmate_coord_optrem.bam --IGNORE=MISSING_TAG_NM --REFERENCE_SEQUENCE=REFGENOME.fasta
```

# Variant calling

[GATK](https://gatk.broadinstitute.org/hc/en-us) version 4.1.9.0

In short: mapping population data to reference - HaplotypeCaller - CombineGVCF - GenotypeGVCFs - SelectVariants

```sh
HaplotypeCaller -R hap0.softmasked.fasta \
        --emit-ref-confidence GVCF \
        -I /T507.sorted.marked_duplicates.bam \
        -O Ppr_gatk_haploT507
```

For each individual
```sh
GATK CombineGVCFs -R hap0.softmasked.fasta \
	-V Ppr_gatk_haploT501 \
	-V Ppr_gatk_haploT502 \
	-V Ppr_gatk_haploT505 \
	-V Ppr_gatk_haploT506 \
	-V Ppr_gatk_haploT507 \
	-O Ppr_merged.g.vcf
```

```sh
GATK GenotypeGVCFs -R hap0.softmasked.fasta \
        -V Ppr_merged.g.vcf \
        -O Ppr_gatk_all.vcf
```

```sh
GATK ValidateVariants --gvcf -R Japan_T504_hifiasm_hap0_scaff10x_softmask.fa \
	-V Ppr_JP_gatk_all.vcf
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
	-O Ppr_gatk_all.snp.f.vcf.gz
```

```sh
GATK VariantsToTable \
	-V Ppr_gatk.SNP.filtered.vcf.gz \
	-F CHROM -F POS -F REF -F ALT -F MULTI-ALLELIC -F TYPE -GF AD \
	-O Ppr_gatk.SNP.filtered.table
```
