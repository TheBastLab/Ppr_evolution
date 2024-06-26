# Genome annotation
	Finally got the annotation file: Ppr.hap0.gff; protein sequenece: Ppr.hap0.pep; CDS file: Ppr.hap0.cds
## software version
	EDTA v1.9.8 	
	Bedtools Version: v2.26.0
	STAR v2.5.1a
	Trinity v2.1.1
	PASA 2.5.2
	stringtie v2.2.0
	TransDecoder 5.5.0
	EVidenceModeler 1.1.1
	BRAKER2 v2.1.6
	gffread v0.12.1
	emapper 2.1.6
	InterProScan version 5.61-93.0
	clusterProfiler v4.6.2

## 1. GenomeMask
### 1.1 Hardmask
	  EDTA.pl --genome polished.fasta --sensitive 1 --anno 1  --threads 50 --overwrite 1
### 1.2 Softmask
	  bedtools maskfasta -fi polished.fasta.masked -fo polished.softmask.fasta -bed polished.fa.mod.EDTA.TEanno.gff3 -soft
	  
## 2. Gene Structure Prediction
### 2.1 RNA-Seq Mapping and Trinity assembly
	STAR \
	--genomeDir STAR \
	--runThreadN 20 \
	--readFilesIn 150360_R1.fq.gz, 150360_R2.fq.gz \
	--readFilesCommand zcat \
	--outFileNamePrefix Ppr_lowmismatch \
	--outSAMtype BAM SortedByCoordinate \
	--outBAMsortingThreadN 10 \
	--outSAMstrandField intronMotif \
	--outFilterIntronMotifs RemoveNoncanonical
 	
	STAR \
	--genomeDir STAR \
	--runThreadN 20 \
	--readFilesIn SRR4039022_1_val_1.fq.gz, SRR4039022_2_val_2.fq.gz \
	--readFilesCommand zcat \
	--outFileNamePrefix Ppr_lowmismatch \
	--outSAMtype BAM SortedByCoordinate \
	--outBAMsortingThreadN 10 \
	--outSAMstrandField intronMotif \
	--outFilterIntronMotifs RemoveNoncanonical

  	Trinity --seqType fq --max_memory 60G --CPU 20  \
  	--left 150360_R1.fq.gz,SRR4039022_1_val_1.fq.gz \
  	--right 150360_R2.fq.gz,SRR4039022_2_val_2.fq.gz \
	--output trinity.out

### 2.2 PASA mapping
	Launch_PASA_pipeline.pl \
	-c alignAssembly.config \
	-C -R -g Ppr.hap0.softmasked.fasta \
	-t Trinity.fasta \
	--ALIGNERS blat --CPU 10

### 2.3 Stringtie assembly
	stringtie -p 50 -o Ppr_hap0_PprAligned_stringtie.gtf PprAligned.sortedByCoord.out.bam
	stringtie -p 50 -o Ppr_hap0_Ppr_downAligned_stringtie.gtf Ppr_downAligned.sortedByCoord.out.bam
	stringtie --merge -p 50  -o merged.gtf Ppr_hap0_PprAligned_stringtie.gtf Ppr_hap0_Ppr_downAligned_stringtie.gtf
	cufflinks_gtf_genome_to_cdna_fasta.pl merged.gtf /RAID/Data/Mites/Genomes/Ppr/German_eiffel/hap0/Ppr.hap0.softmasked.fasta > Ppr_merged.fasta
	cufflinks_gtf_to_alignment_gff3.pl merged.gtf > transcripts.gff3
	TransDecoder.LongOrfs -t Ppr_merged.fasta > TransDecoder.LongOrfs.log
	TransDecoder.Predict -t Ppr_merged.fasta > TransDecoder.Predict.log 
	cdna_alignment_orf_to_genome_orf.pl Ppr_merged.fasta.transdecoder.gff3 transcripts.gff3 Ppr_merged.fasta > Ppr_merged.fasta.transdecoder.genome.gff3

### 2.4 BRAKER2(AUGUEST)
	braker.pl --cores 40 --species=Ppr_2bam --genome=Ppr.hap0.softmasked.fasta \
		--softmasking --bam=PprAligned.sortedByCoord.out.bam, Ppr_downAligned.sortedByCoord.out.bam \
		--gff3 --useexisting --workingdir=braker_out
    
### 2.5 EVidenceModeler
	ln /pasa.sqlite.pasa_assemblies.gff3 transcript_alignments.gff3
	cat Ppr_merged.fasta.transdecoder.genome.gff3 augustus.hints.gff3 > gene_predictions.gff3
	export genome=Ppr.hap0.softmasked.fasta
	export prefix=Ppr_hap0
	partition_EVM_inputs.pl --genome $genome --gene_predictions `pwd`/gene_predictions.gff3 --transcript_alignments `pwd`/transcript_alignments.gff3 --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out
  
  	write_EVM_commands.pl --genome $genome --weights weights.txt \
	  --gene_predictions `pwd`/gene_predictions.gff3 \
	  --transcript_alignments `pwd`/transcript_alignments.gff3 \
	  --output_file_name evm.out  --partitions partitions_list.out >  commands.list
	parallel --jobs 10 < commands.list
  
	recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out
	convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome $genome
	find . -regex ".*evm.out.gff3" -exec cat {} \; | bedtools sort -i - > EVM.all.gff
  
### 2.6 Protein and coding sequence 
	gffread $prefix.gff -g $genome -y $prefix.$prefix.gff.pep
	gffread $prefix.gff -g $genome -x $prefix.$prefix.gff.cds

## 3. Functional Annotation
### 3.1 EggNOG-mapper 
	emapper.py --cpu 10 -i hap0_UTR_longest.pep -o hap0_UTR_longest
 
### 3.2 InterProScan
	interproscan.sh -i hap0_UTR_longest.pep -T ./tmp -goterms -iprlookup -pa

### 3.3 Merge Go annotation
	python GO_anno_from_tab.py -i hap0_UTR_longest.pep_1.tsv
	python GO_anno_from_tab.py -i hap0_UTR_longest.emapper.annotations

### 3.4 clusterProfiler(working in R)
	require(clusterProfiler)
	egg<-read.delim("merged_Go")
	data <- read.table("hap1tohap2.diff_rpkm1,header=F)
	genes <- as.character(data$V1)
	dea <- enricher(genes, TERM2GENE = egg[, c("GO_term", "Gene_ID")], TERM2NAME = egg[, c("GO_term", "Function")], pAdjustMethod= "BH",pvalueCutoff  = 0.05, qvalueCutoff  = 0.05)
	write.table(as.data.frame(dea),"go_dea_enrich.csv",sep="\t",row.names =F,quote=F)



