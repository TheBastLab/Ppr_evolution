## software
gatk=/NVME/Software/popgen/gatk-4.1.9.0/gatk

## data
ref=/home/hoeztopr/Data/hoeztopr/Ppr/genomes/DE/Ppr.hap0.softmasked.fasta
bam=Ppr_T501.sorted.bam #generated with mapping.sh


## Preparation
#### build dict for genome
```
/NVME/Software/popgen/gatk-4.1.9.0/gatk CreateSequenceDictionary \
-R /home/hoeztopr/Data/hoeztopr/Ppr/genomes/DE/Ppr.hap0.softmasked.fasta -O Ppr.dict
```
```
samtools faidx Ppr.russia.hap0.softmasked.fasta
```
#### generate bam files 
```
/home/hoeztopr/Data/hoeztopr/Scripts/mapping.sh
```
        #### deduplication removing with picard
        /home/hoeztopr/Data/hoeztopr/Scripts/deduplication.sh
        samtools index -@ 60 T506.sorted.marked_duplicates.bam

        ##### compare duplication removed bam with regular bam
        /home/hoeztopr/Data/hoeztopr/Scripts/comparebams.sh
### Check bamfiles
```
/NVME/Software/popgen/gatk-4.1.9.0/gatk GetSampleName \
     -I /RAID/Data/gaoshan/gaoshan/hifiasm_tell-seq/tell_sort/tell_sort_out/Ppr_T502/Ppr_T502_temp/Ppr_T502.sorted.bam \
     -O sample_nameT502.txt
```
### check if bam is good to go for GATK ###
```
gatk ValidateSamFile --INPUT=INDIVIDUAL_readgroup_fixmate_coord_optrem.bam --IGNORE=MISSING_TAG_NM --REFERENCE_SEQUENCE=REFGENOME.fasta
```
### Prep the Genome
```
/NVME/Software/popgen/gatk-4.1.9.0/gatk CreateSequenceDictionary -R Ppr_instagrall.polished.FINAL.softmask.fa
```

# Call gvcf per sample
### DE pop: T501,2,5,6,7
```
/NVME/Software/popgen/gatk-4.1.9.0/gatk HaplotypeCaller \
        -R /home/hoeztopr/Data/hoeztopr/Ppr/genomes/DE/Ppr.hap0.softmasked.fasta \
        --emit-ref-confidence GVCF \
        -I /home/hoeztopr/Data/hoeztopr/Ppr/mapped/DE/T507.sorted.marked_duplicates.bam \
        -O /home/hoeztopr/Data/hoeztopr/Ppr/gvcfs/DE/Ppr_gatk_haploT507
                  #####  --> for all individuals
                  ##### --native-pair-hmm-threads 4 (8-10) OR run simultaneasly
```
# Merge them
```
/NVME/Software/popgen/gatk-4.1.9.0/gatk CombineGVCFs \
-R /home/hoeztopr/Data/hoeztopr/Ppr/genomes/DE/Ppr.hap0.softmasked.fasta \
-V Ppr_gatk_haploT501 \
-V Ppr_gatk_haploT502 \
-V Ppr_gatk_haploT505 \
-V Ppr_gatk_haploT506 \
-V Ppr_gatk_haploT507 \
-O Ppr_merged.g.vcf
```
##### Processed 101.962.705 total variants in 27.8 minutes.

# Detect SNPs
```
/NVME/Software/popgen/gatk-4.1.9.0/gatk GenotypeGVCFs \
        -R /home/hoeztopr/Data/hoeztopr/Ppr/genomes/DE/Ppr.hap0.softmasked.fasta \
        -V Ppr_merged.g.vcf \
        -O Ppr_gatk_all.vcf
```
##### 7.202.661 total variants in 30.6 minutes.

        # Validate variants
        ```
        /home/hoeztopr/Software/popgen/gatk-4.2.2.0/gatk ValidateVariants \
        --gvcf \
        -R /home/hoeztopr/Data/hoeztopr/Ppr/genomes/DE/Ppr.hap0.softmasked.fasta \
        -V Ppr_gatk_all.vcf
        ```
# Compress files
```
bgzip -f Ppr_gatk_all.vcf
tabix -p vcf Ppr_gatk_all.vcf.gz
```
```
rm Ppr_merged.g.vcf
```
# Filter out only SNPs from VCF
```
/NVME/Software/popgen/gatk-4.1.9.0/gatk SelectVariants \
-select-type SNP \
-V Ppr_gatk_all.vcf.gz \
-O Ppr_gatk_all.snp.vcf.gz
```
##### Processed 7.038.890 total variants

# Hardfiltering SNPs by parameters
```
/NVME/Software/popgen/gatk-4.1.9.0/gatk VariantFiltration \
-V Ppr_gatk_all.snp.vcf.gz  \
--filter-expression "QD < 2.0" --filter-name "QD2" \
--filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
--filter-expression "SOR > 3.0" --filter-name "SOR3" \
--filter-expression "FS > 60.0" --filter-name "FS60" \
--filter-expression "MQ < 40.0" --filter-name "MQ40" \
--filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
--output Ppr_gatk.SNP.filtered.vcf.gz
```
##### Processed 6745225 total variants in 4.0 minutes.





--------------------------------------------------------------------
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
