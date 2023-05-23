# Differential expression of alleles
## software version
    STAR v2.5.1a
    GFOLD V1.1.4
    samtools v1.11
## 1. RNA-Seq Mapping
    STAR \
    --genomeDir STAR \
    --runThreadN 40 \
    --readFilesIn 150360_R1.fq.gz, 150360_R2.fq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix Ppr_lowmismatch \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMsortingThreadN 10 \
    --outSAMstrandField intronMotif \
    --outFilterMismatchNmax 3 \
    --outFilterMismatchNoverLmax 0.1 \
    --outFilterMismatchNoverReadLmax 0.5 \
    --outFilterIntronMotifs RemoveNoncanonical

## 2. GFOLD calculate different expression genes
    samtools view PprAligned.sortedByCoord.out.bam | gfold count -ann hap1_add_UTR.gtf3 -tag stdin -o hap1.read_cnt
    samtools view PprAligned.sortedByCoord.out.bam | gfold count -ann hap2_add_UTR.gtf3 -tag stdin -o hap2.read_cnt
    gfold diff -s1 hap1 -s2 hap2 -suf .read_cnt -o hap1VShap2.diff
    
## 3. Filter GFOLD and RPKM value
    awk '($7+$6>1)&&($3>0.05||$3<-0.05){print$0}' hap1tohap2.diff > hap1tohap2.diff_rpkm1_gfold0.05.genelist
