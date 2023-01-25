# *Platynothrus peltifer* genome analysis
*Platynothrus peltifer* genome with haplotype-specific analyses 

## Table of contents
* [Genome assembly pipeline](#Genome-assembly-pipeline)
* [Genome annotation](#Genome-annotation)
* [Variant calling](#Variant-calling)

## Genome assembly pipeline

### *k*-mer analysis

[KAT](https://github.com/TGAC/KAT) version 2.4.2
```sh
kat hist -o kat_hist hifi_reads.fastq.gz
kat gcp -o kat_gcp hifi_reads.fastq.gz
```

<img src="./fig/kat_hist.svg" width=400> <img src="./fig/kat_gcp.svg" width=400>

[Smudgeplot](https://github.com/KamilSJaron/smudgeplot) version 0.2.5
[KMC](https://github.com/tbenavi1/KMC) version 
```sh
mkdir tmp_smudge
ls hifi_reads.fastq.gz > FILES
kmc -k27 -ci1 -cs10000 @FILES kmcdb tmp_smudge
kmc_tools transform kmcdb histogram kmcdb_k27.hist -cx10000

L=$(smudgeplot.py cutoff kmcdb_k27.hist L)
U=$(smudgeplot.py cutoff kmcdb_k27.hist U)
echo $L $U

kmc_tools transform kmcdb -ci"$L" -cx"$U" dump -s kmcdb_L"$L"_U"$U".dump
smudgeplot.py hetkmers -o kmcdb_L"$L"_U"$U" < kmcdb_L"$L"_U"$U".dump

smudgeplot.py plot kmcdb_L"$L"_U"$U"_coverages.tsv
```

### *De novo* assembly

[hifiasm](https://github.com/chhylp123/hifiasm) version 0.16.1-r375
```sh
hifiasm -o assembly hifi_reads.fastq.gz
```

### TELL-seq linked reads improve assembly by Scaff10X
	/home/shangao/Software/Assembly/Scaff10X/src/scaff10x \
    -nodes 30 -longread 1 -gap 100 -matrix 2000 -reads 10 -score 10 -edge 50000 -link 8 -block 50000 -plot Ppr_result/barcode_lengtg.png \
    /home/shangao/Scratch/gaoshan/hifiasm/genome/Ppr_new/Ppr_new.fa \
    genome-BC_1.fastq.gz \
    genome-BC_2.fastq.gz \
    output_scaffolds.fastae

### Omni-C scaffolding 

[bwa](https://github.com/lh3/bwa) version 0.7.15
[hicstuff](https://github.com/koszullab/hicstuff) version 3.1.1

```sh
hicstuff pipeline -e 100 -a bwa -g assembly.fasta -m iterative \
	-o hicstuff_out omnic.trimmed.end1.fastq omnic.trimmed.end2.fastq
```

[instaGRAAL](https://github.com/koszullab/instaGRAAL) version 0.1.6 no-opengl branch

```sh
instagraal --level 5 --cycles 100 hicstuff_out assembly.fasta instagraal_out
```

```sh
instagraal-polish -m polishing -f assembly.purged.fasta -j NNNNNNNNNN \
	-i instagraal_out/hicstuff_out/test_mcmc5/info_frags.txt \
	-o assembly.hic_scaffolds.fasta
```
	
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

### Assembly evaluation


[KAT](https://github.com/TGAC/KAT) version 2.4.2
```sh
kat comp -o kat_comp_hap0 hifi_reads.fastq.gz hap0.fasta
cat hapA.fasta hapB.fasta > phased.fasta
kat comp -o kat_comp_phased hifi_reads.fastq.gz phased.fasta
kat comp -o kat_comp_hapA hifi_reads.fastq.gz hapA.fasta
kat comp -o kat_comp_hapB hifi_reads.fastq.gz hapB.fasta
```


<img src="./fig/figure_kmer_comp.svg" width=800>

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

1. GATK: /home/hoeztopr/Data/hoeztopr/Scripts/vcfcall.txt
software: /NVME/Software/popgen/gatk-4.1.9.0/gatk
- In short: mapping population data to reference - HaplotypeCaller - CombineGVCF - GenotypeGVCFs - SelectVariants

```sh
HaplotypeCaller -R /home/hoeztopr/Data/hoeztopr/Ppr/genomes/DE/Ppr.hap0.softmasked.fasta \
        --emit-ref-confidence GVCF \
        -I /home/hoeztopr/Data/hoeztopr/Ppr/mapped/DE/T507.sorted.marked_duplicates.bam \
        -O /home/hoeztopr/Data/hoeztopr/Ppr/gvcfs/DE/Ppr_gatk_haploT507
```
for each individual
```sh
/NVME/Software/popgen/gatk-4.1.9.0/gatk CombineGVCFs \
-R /home/hoeztopr/Data/hoeztopr/Ppr/genomes/DE/Ppr.hap0.softmasked.fasta \
-V Ppr_gatk_haploT501 \
-V Ppr_gatk_haploT502 \
-V Ppr_gatk_haploT505 \
-V Ppr_gatk_haploT506 \
-V Ppr_gatk_haploT507 \
-O Ppr_merged.g.vcf
```
```sh
/NVME/Software/popgen/gatk-4.1.9.0/gatk GenotypeGVCFs \
        -R /home/hoeztopr/Data/hoeztopr/Ppr/genomes/DE/Ppr.hap0.softmasked.fasta \
        -V Ppr_merged.g.vcf \
        -O /home/hoeztopr/Data/hoeztopr/Ppr/vcf/Ppr_gatk_all.vcf
```
```sh
/home/hoeztopr/Software/popgen/gatk-4.2.2.0/gatk ValidateVariants \
--gvcf \
-R /home/hoeztopr/Data/hoeztopr/Ppr/genomes/Japan_T504_hifiasm_hap0_scaff10x_softmask.fa \
-V Ppr_JP_gatk_all.vcf
```
```sh
/NVME/Software/popgen/gatk-4.1.9.0/gatk SelectVariants \
-select-type SNP \
-V Ppr_gatk_all.vcf.gz \
-O Ppr_gatk_all.snp.vcf.gz
```
```sh
/NVME/Software/popgen/gatk-4.1.9.0/gatk VariantFiltration \
-V Ppr_gatk_all.snp.vcf.gz \
--filter-expression "QD <2.0 || MQ <40.0 || FS >60.0 || SOR >3.0 || ReadPosRankSum < -8.0 || MQRankSum < -12.5" \
--filter-name "PASS" \
-O Ppr_gatk_all.snp.f.vcf.gz
```
```sh
/NVME/Software/popgen/gatk-4.1.9.0/gatk VariantsToTable \
-V Ppr_gatk.SNP.filtered.vcf.gz \
-F CHROM -F POS -F REF -F ALT -F MULTI-ALLELIC -F TYPE -GF AD \
-O Ppr_gatk.SNP.filtered.table
```


