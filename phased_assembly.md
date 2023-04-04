# Phased assembly pipeline 

## *De novo* assembly

[Flye](https://github.com/fenderglass/Flye) version 2.9
```sh
flye -o flye_v29_default --pacbio-raw hifi_reads.fastq.gz 
```

## TELL-Seq scaffolding

Tell-Read version 1.0.3
```
mkdir -p data/genome/phased
cd data/genome/phased
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
scaff10x -longread 1 -gap 100 -matrix 2000 -reads 10 -score 10 -edge 50000 -link 8 -block 50000 -plot barcode_length.png draft_assembly.fasta tellseq.end1.fastq.gz tellseq.end2.fastq.gz scaff10x_scaffolds.fasta
```

## Haplotype separation

[minimap2](https://github.com/lh3/minimap2) version 2.24-r1122
[purge_dups](https://github.com/dfguan/purge_dups) version 1.2.5
```sh
pri_asm="scaffolded"

minimap2 -x map-hifi ${pri_asm}.fasta hifi_reads.fastq.gz | gzip -c - > minimap2_hifi.scaffolded.paf.gz

split_fa $pri_asm > $pri_asm.split
minimap2 -xasm5 -DP $pri_asm.split $pri_asm.split | gzip -c - > $pri_asm.split.self.paf.gz

pbcstat minimap2_hifi.${pri_asm}.paf.gz 
calcuts -l 5 -m 73 -u 327 PB.stat > cutoffs 2>calcults.log

hist_plot.py -c cutoffs PB.stat purge_dups.${pri_asm}.png

purge_dups -2 -T cutoffs -c PB.base.cov $pri_asm.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs dups.bed $pri_asm
mv purged.fa scaffolded.alt1.fasta
mv hap.fa scaffolded.alt2.fasta
```

## RagTag scaffolding

[RagTag](https://github.com/malonge/RagTag) version 2.1.0
```sh
ragtag.py scaffold -o ragtag_alt1_scaffolded_alt2 scaffolded.alt2.fasta scaffolded.alt1.fasta
ragtag.py scaffold -o ragtag_alt2_scaffolded_alt2 scaffolded.alt1.fasta scaffolded.alt2.fasta
```
