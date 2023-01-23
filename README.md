# *Platynothrus peltifer* genome analysis
*Platynothrus peltifer* genome with haplotype-specific analyses 

## Table of contents
* [Genome assembly pipeline](#Genome-assembly-pipeline)

## 1. Genome assembly pipeline

### 1.1 *k*-mer analysis

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

### 1.2 *De novo* assembly

[hifiasm](https://github.com/chhylp123/hifiasm) version 0.16.1-r375
```sh
/home/shangao/software/hifiasm-0.16.1/hifiasm -o test -t 50 /home/shangao/Data/../mites/reads/PacBio/Ppr/m64093_200831_134054.Q20.fastq.gz
```

### 1.3 Tell-Seq linked reads improve assembly by Scaff10X
	/home/shangao/Software/Assembly/Scaff10X/src/scaff10x \
    -nodes 30 -longread 1 -gap 100 -matrix 2000 -reads 10 -score 10 -edge 50000 -link 8 -block 50000 -plot Ppr_result/barcode_lengtg.png \
    /home/shangao/Scratch/gaoshan/hifiasm/genome/Ppr_new/Ppr_new.fa \
    genome-BC_1.fastq.gz \
    genome-BC_2.fastq.gz \
    output_scaffolds.fastae

### 1.4 Omni-C scaffolding 

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

### 1.5 Genome quality check
#### 1.5.1 N50
	~/Software/tools/bbmap/bbstats.sh in=output_scaffolds.fasta > bbstats.out
	
#### 1.5.2 Busco
	busco -i .fa -m genome -o busco --auto-lineage-euk
	
#### 1.5.3 Blobtools2
	conda activate btk_env
	
##### Create database
	/home/shangao/Software/blobtoolkit/blobtools2/blobtools create \
    --fasta .fa \
    /home/shangao/Scratch/gaoshan/hifiasm/test_Ppr/blobtools2/dataset/Ppr
    
##### Add coverage, hits, busco results
	/home/shangao/Software/blobtoolkit/blobtools2/blobtools add \
	--taxrule bestsumorder \
	--taxdump /RAID/Data/databases/taxdump/ \
	--cov /home/shangao/Scratch/gaoshan/hifiasm/test_Ppr/Ppr_Pacbio1.bam \
	--hits /home/shangao/Scratch/gaoshan/hifiasm/test_Ppr/Ppr.ncbi.blastn.out \
  	--busco /home/shangao/Scratch/hifiasm/test_Ppr/busco_p/run_arachnida_odb10/full_table.tsv \
    	--busco /home/shangao/Scratch/hifiasm/test_Ppr/busco_p/run_eukaryota_odb10/full_table.tsv \
	/home/shangao/Scratch/gaoshan/hifiasm/test_Ppr/blobtools2/dataset/Ppr

##### Filter
	/home/shangao/Software/blobtoolkit/blobtools2/blobtools filter \
    	--param length--Min=1000 \
	--summary-rank species \
 	--fasta ../test.p_ctg.fa \
	/home/shangao/Scratch/gaoshan/hifiasm/test_Ppr/blobtools2/dataset/Ppr
##### View
	/home/shangao/Software/blobtoolkit/blobtools2/blobtools view --remote /home/shangao/Scratch/hifiasm/test_Ppr/blobtools2/dataset/Ppr
	
### 1.6 Gap filling
Version : 1.1.1

	/home/shangao/Software/TGS-GapCloser/TGS-GapCloser.sh \
	--scaff ../Ppr_instagrall.FINAL.sort.fasta \
	--reads /home/shangao/Data/../mites/reads/PacBio/Ppr/m64093_200831_134054.Q20.fastq.gz \
	--output Ppr_without_correct \
	--tgstype pb \
        --ne \
	--minmap_arg '-x asm20'

### 1.7 Polishing 
v1.0.3

	hypo -d Ppr_without_correct.scaff_seqs \
	-r /home/shangao/Data/../mites/reads/PacBio/Ppr/m64093_200831_134054.Q20.fastq.gz \
	-s 200m -c 100 -b mapped-ccs.sorted.bam -t 48 \
	-o Ppr_instagrall.polished.fa

## 2. Genome annotation
	








## TO DO
1. please add all paths to the basic data

## VCF calling
Data
- reference: Ppr_instagrall.polished.FINAL.softmask.fa
- bamfiles: /home/hoeztopr/Data/hoeztopr/VCFcalling/GATK_version/bamfilenames.list (/RAID/Data/gaoshan/hifiasm_tell-sort/)
- haplotypes: /home/hoeztopr/Data/hoeztopr/VCFcalling/GATK_version #generated with GATK HaplotypeCaller

1. GATK: /home/hoeztopr/Data/hoeztopr/Scripts/vcfcall.txt
software: /NVME/Software/popgen/gatk-4.1.9.0/gatk
- In short: CombineGVCF - GenotypeGVCFs - SelectVariants

Path_to_data: /home/hoeztopr/Data/hoeztopr/VCFcalling/GATK_version

2. Freebayes: /home/hoeztopr/Data/hoeztopr/Scripts/Freebayes.txt
software= ~/Software/popgen/freebayes1.3.0/freebayes-v1.3.0-1
- In short: version a) simple call variants -> Ppr.fb.vcf b) with population information -> Ppr.pop.fb.vcf 
- Filtering: GATK SelectVariants(Ppr.fb.gatk.snp.gz) + Hardfiltering:VariantFiltration (Ppr.fb.gatk.snp.f.gz)

Path_to_data: /home/hoeztopr/Data/hoeztopr/VCFcalling/Freebayes

3. TELL-sort
VCF provided by Shan: Ppr.merged.vcf.gz
- Filtering: VCFtools (Ppr.SHAN.filtered.recode.vcf) # Hardfiltering: Ppr.SHAN.filtered_GATKhard.vcf.gz

Path_to_data: /home/hoeztopr/Scratch/hoeztopr/nucdiv/SHAN_VCF
123
