# Assembly evaluation

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
