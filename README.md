# *Platynothrus peltifer* genome analysis
*Platynothrus peltifer* genome with haplotype-specific analyses 

## Table of contents
* [Initial analyses on HiFi reads](#Initial-analyses-on-HiFi-reads)
* [Collapsed chromosome-level assembly pipeline for German Ppr](https://github.com/TheBastLab/Ppr_evolution/blob/main/assembly_hap0_German_Ppr.md)
* [Collapsed assembly pipeline for Italian and Russian Ppr](https://github.com/TheBastLab/Ppr_evolution/blob/main/initial_hifi_analyses.md)
* [Phased assembly pipeline](#Phased-assembly-pipeline)
* [Assembly evaluation](#Assembly-evaluation)
* [Genome annotation](#Genome-annotation)
* [Variant calling](#Variant-calling)



## Collapsed assembly pipeline for Italian and Russian Ppr

### *De novo* assembly

[Flye]() v2.9
```sh
flye -o flye_v29_default --pacbio-raw hifi_reads.fastq.gz 
```

[hifiasm](https://github.com/chhylp123/hifiasm) version 0.16.1-r375
```sh
hifiasm -o assembly hifi_reads.fastq.gz
```

### Haplotig purging

[purge_dups](https://github.com/dfguan/purge_dups) v1.2.5
```sh

```

### TELL-seq scaffolding 

### RagTag scaffolding

### Gap filling

[TGS-GapCloser](https://github.com/BGI-Qingdao/TGS-GapCloser) version 1.1.1

```sh
TGS-GapCloser.sh --scaff scaffolds.fasta \
	--reads hifi_reads.fastq.gz \
	--output gap_filled \
	--tgstype pb --ne \
	--minmap_arg '-x asm20'
```

### Polishing 

[HyPo](https://github.com/kensung-lab/hypo) v1.0.3
[minimap2](https://github.com/lh3/minimap2) version 2.24r1122
[SAMtools](https://github.com/samtools/samtools) version 1.11
```sh
minimap2 --secondary=no --MD -ax map-hifi gap_filled.fasta hifi_reads.fastq.gz | samtools view -Sb - > mapped-ccs.bam
samtools sort -o mapped-ccs.sorted.bam mapped-ccs.bam
samtools index mapped-ccs.sorted.bam

hypo -d gap_filled.fasta -r hifi_reads.fastq.gz -s 200m -c 100 -b mapped-ccs.sorted.bam \
	-o polished.fasta
```

## Phased assembly pipeline 

### *De novo* assembly

[Flye]() v2.9
```sh
flye -o flye_v29_default --pacbio-raw hifi_reads.fastq.gz 
```

### TELL-Seq scaffolding

### Haplotype separation

[purge_dups](https://github.com/dfguan/purge_dups) v1.2.5
```sh

```

### RagTag scaffolding



## Assembly evaluation

[KAT](https://github.com/TGAC/KAT) version 2.4.2
```sh
kat comp -o kat_comp_hap0 hifi_reads.fastq.gz hap0.fasta
cat hapA.fasta hapB.fasta > phased.fasta
kat comp -o kat_comp_phased hifi_reads.fastq.gz phased.fasta
kat comp -o kat_comp_hapA hifi_reads.fastq.gz hapA.fasta
kat comp -o kat_comp_hapB hifi_reads.fastq.gz hapB.fasta
```

[BUSCO](https://busco.ezlab.org/) version 5.0.0

```sh
busco -i final_scaffolds.fasta -m genome -o busco_out_arachnida_odb10 -l arachnida_odb10
busco -i final_scaffolds.fasta -m genome -o busco_out_arthropoda_odb10 -l arthropoda_odb10
```

[minimap2](https://github.com/lh3/minimap2) version 2.24r1122
[SAMtools](https://github.com/samtools/samtools) version 1.11
```sh
minimap2 -ax map-hifi hap0.fasta hifi_reads.fastq.gz | samtools view -b | samtools sort -o minimap2_hifi.hap0.bam
```

[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) version 2.6.0
```sh
blastn -query hap0.fasta -db nt -outfmt "6 qseqid staxids bitscore std sscinames scomnames" \
	-max_hsps 1 -evalue 1e-25 -out blast.out
```

[BlobTools2](https://blobtoolkit.genomehubs.org/blobtools2/) version 2.3.3
```sh
blobtools add --fasta hap0.fasta --cov minimap2_hifi.hap0.bam --hits blast.out \
	--busco busco_out_arachnida_odb10/run_arachnida_odb10/full_table.tsv \
	--taxdump taxdump --create hap0_BLOBDIR
blobtools view hap0_BLOBDIR
```

<img src="./fig/final_scaffolds_BLOBDIR.blob.circle.svg" width=600>

## Genome annotation

## Variant calling

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


