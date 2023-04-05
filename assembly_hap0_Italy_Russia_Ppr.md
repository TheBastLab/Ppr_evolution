
# Collapsed assembly pipeline for Italian and Russian Ppr

## *De novo* assembly

[Flye]() v2.9
```sh
flye -o flye_v29_default --pacbio-raw hifi_reads.fastq.gz 
```

[hifiasm](https://github.com/chhylp123/hifiasm) version 0.16.1-r375
```sh
hifiasm -o assembly hifi_reads.fastq.gz
```

## Haplotig purging

For both flye and hifiasm assemblies.

[minimap2](https://github.com/lh3/minimap2) version 2.24-r1122
[purge_dups](https://github.com/dfguan/purge_dups) version 1.2.5
```sh
pri_asm="draf"

minimap2 -x map-hifi ${pri_asm}.fasta hifi_reads.fastq.gz | gzip -c - > minimap2_hifi.scaffolded.paf.gz

split_fa $pri_asm > $pri_asm.split
minimap2 -xasm5 -DP $pri_asm.split $pri_asm.split | gzip -c - > $pri_asm.split.self.paf.gz

pbcstat minimap2_hifi.${pri_asm}.paf.gz 
calcuts -l L -m M -u U PB.stat > cutoffs 2>calcults.log

hist_plot.py -c cutoffs PB.stat purge_dups.${pri_asm}.png

purge_dups -2 -T cutoffs -c PB.base.cov $pri_asm.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs dups.bed $pri_asm
mv purged.fa scaffolded.alt1.fasta
mv hap.fa scaffolded.alt2.fasta
```

## TELL-seq scaffolding 

For both Flye+purge_dups and hifiasm+purge_dups assemblies.

Tell-Read version 1.0.3
```
mkdir -p data/genome/hap0
cd data/genome/hapO
cp draft_assembly.fasta .
generateGenomeIndexBed.sh draft_assembly.fasta

cd .. 
mv ../genomes.json .
mkdir data/TELL-seq_preprocessing

run_tellread.sh -i sequencing_data/ -o data/TELL-seq_preprocessing -f data/genome -s BARCODE -g LABEL

ust10x -sz 200000000 \
	-i1 data/TELL-seq_preprocessing/1_demult/Full/TELL-seq_preprocessing_I1.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz \
	-r1 data/TELL-seq_preprocessing/1_demult/Full/TELL-seq_preprocessing_R1.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz \
	-r2 data/TELL-seq_preprocessing/1_demult/Full/TELL-seq_preprocessing_R2.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz \
	-wl 4M-with-alts-february-2016.txt 
```

[Scaff10X](https://github.com/wtsi-hpag/Scaff10X) version 4.2
```
scaff_reads ./input.dat tellseq.end1.fastq.gz tellseq.end2.fastq.gz
scaff10x -longread 1 -gap 100 -matrix 2000 -reads 10 -score 10 -edge 50000 -link 8 -block 50000 -plot barcode_length.png draft_assembly.fasta tellseq.end1.fastq.gz tellseq.end2.fastq.gz tellseq_scaffolds.fasta
```

## RagTag scaffolding

[RagTag](https://github.com/malonge/RagTag) version 2.1.0
```sh
ragtag.py scaffold -o ragtag_flye_scaffolded_hifiasm hifiasm.purge_dups.fatsa flye.purge_dups.fasta 
```

## Gap filling

[TGS-GapCloser](https://github.com/BGI-Qingdao/TGS-GapCloser) version 1.1.1

```sh
TGS-GapCloser.sh --scaff scaffolds.fasta --reads hifi_reads.fastq.gz --output gap_filled --tgstype pb --ne --minmap_arg '-x asm20'
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
