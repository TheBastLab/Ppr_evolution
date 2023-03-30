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
