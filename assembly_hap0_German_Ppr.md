# Collapsed chromosome-level assembly pipeline

## *De novo* assembly

[hifiasm](https://github.com/chhylp123/hifiasm) version 0.16.1-r375
```sh
hifiasm -o assembly hifi_reads.fastq.gz
```

## Linked-read scaffolding

[Scaff10X](https://github.com/wtsi-hpag/Scaff10X) version 4.2
```sh
scaff10x -longread 1 -gap 100 -matrix 2000 -reads 10 -score 10 -edge 50000 -link 8 -block 50000 -plot barcode_length.png \
    assembly.fasta genome-BC_1.fastq.gz genome-BC_2.fastq.gz  linked_scaffolds.fasta
```

## Omni-C scaffolding 

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
	
## Gap filling

[TGS-GapCloser](https://github.com/BGI-Qingdao/TGS-GapCloser) version 1.1.1

```sh
TGS-GapCloser.sh --scaff scaffolds.fasta \
	--reads hifi_reads.fastq.gz \
	--output gap_filled \
	--tgstype pb --ne \
	--minmap_arg '-x asm20'
```

## Polishing 

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
