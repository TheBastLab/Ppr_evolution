#competitive mapping
for i in $(cat list); do bbsplit.sh threads=40 ref=hap1/Eiffel_hap1.fa,hap2/Eiffel_hap2.fa in=eifrusjap/$i\_R1.fastq.gz in2=eifrusjap/$i\_R2.fastq.gz basename=bbsplit_split/$i-%-bbsplit_R#.fastq.gz refstats=bbsplit_split/$i-refstats scafstats=bbsplit_split/$i-scafstats ambiguous2=split minratio=0.9; done


#mapping of split sets of resequencing reads (shown for HapA only)
bwa index ../hap0/Eiffel_hap0.fa
for i in $(cat ../list); do \
bwa mem -t 40 -R "@RG\tID:${i}_1\tSM:${i}_1\tLB:${i}_1\tPL:Illumina" \
../hap0/Eiffel_hap0.fa \
${i}-Eiffel_hap1-bbsplit_R1_all.fastq.gz \
${i}-Eiffel_hap1-bbsplit_R2_all.fastq.gz > ${i}_alt1.sam
samtools view -bS ${i}_alt1.sam > ${i}_alt1.bam
samtools sort ${i}_alt1.bam -o ${i}_alt1.sorted.bam
samtools index -@ 40 ${i}_alt1.sorted.bam 
rm -f ${i}_alt1.sam
rm -f ${i}_alt1.bam
done


#deduplication (shown for HapA only)
for i in $(cat ../list); do \
picard MarkDuplicates \
I=${i}_alt1.sorted.bam \
O=${i}_alt1.sorted.marked_duplicates.bam \
M=${i}_alt1.marked_dup_metrics.txt \
REMOVE_DUPLICATES=true \
REMOVE_SEQUENCING_DUPLICATES=true
done


#variant calling for each deduplicated alignment of split sequencing reads (each haplotype)
gatk CreateSequenceDictionary -R ../hap0/Eiffel_hap0.fa
samtools faidx ../hap0/Eiffel_hap0.fa
for i in $(cat ../list); do samtools index ${i}_alt1.sorted.marked_duplicates.bam; done
gatk HaplotypeCaller -I Eiffel_T501_alt1.sorted.marked_duplicates.bam -O gvcf/Eiffel_T501_1.gvcf.gz -R ../hap0/Eiffel_hap0.fa -ERC BP_RESOLUTION --sample-name Eiffel_T501_1 --sample-ploidy 2 --max-reads-per-alignment-start 0

gatk CombineGVCFs \
	--java-options "-Xmx1000G" \
	-V G_01_A.gvcf.gz \
	-V G_02_A.gvcf.gz \
        -V G_05_A.gvcf.gz \
        -V G_06_A.gvcf.gz \
        -V G_07_A.gvcf.gz \
        -V G_01_B.gvcf.gz \
        -V G_02_B.gvcf.gz \
        -V G_05_B.gvcf.gz \
        -V G_06_B.gvcf.gz \
        -V G_07_B.gvcf.gz \
        -V J_01_A.gvcf.gz \
        -V J_02_A.gvcf.gz \
        -V J_03_A.gvcf.gz \
        -V J_04_A.gvcf.gz \
        -V J_05_A.gvcf.gz \
        -V J_01_B.gvcf.gz \
        -V J_02_B.gvcf.gz \
        -V J_03_B.gvcf.gz \
        -V J_04_B.gvcf.gz \
        -V J_05_B.gvcf.gz \
        -V R_06_A.gvcf.gz \
        -V R_07_A.gvcf.gz \
        -V R_09_A.gvcf.gz \
        -V R_10_A.gvcf.gz \
        -V R_11_A.gvcf.gz \
        -V R_06_B.gvcf.gz \
        -V R_07_B.gvcf.gz \
        -V R_09_B.gvcf.gz \
        -V R_10_B.gvcf.gz \
        -V R_11_B.gvcf.gz \
  	-V I_01_A.gvcf.gz \
        -V I_02_A.gvcf.gz \
        -V I_22_A.gvcf.gz \
        -V I_23_A.gvcf.gz \
        -V I_24_A.gvcf.gz \
        -V I_01_B.gvcf.gz \
        -V I_02_B.gvcf.gz \
        -V I_22_B.gvcf.gz \
        -V I_23_B.gvcf.gz \
        -V I_24_B.gvcf.gz \
        -V C_04_A.gvcf.gz \
        -V C_05_A.gvcf.gz \
        -V C_06_A.gvcf.gz \
        -V C_07_A.gvcf.gz \
        -V C_08_A.gvcf.gz \
        -V C_04_B.gvcf.gz \
        -V C_05_B.gvcf.gz \
        -V C_06_B.gvcf.gz \
        -V C_07_B.gvcf.gz \
        -V C_08_B.gvcf.gz \
	-R ../../hap0/Eiffel_hap0.fa \
	-O merge_AB.gvcf.gz

gatk GenotypeGVCFs --include-non-variant-sites -R ../../hap0/Eiffel_hap0.fa -V merge_AB.gvcf.gz -O geno_AB.gvcf.gz


#filter variants
bcftools view --exclude-types indels --threads 40 -O z geno_AB.gvcf.gz > geno_AB_snp.gvcf.gz
tabix -p gzvcf geno_AB_snp.gvcf.gz
bcftools view -m 1 -M 2 --threads 40 -O z geno_AB_snp.gvcf.gz > geno_AB_snp_an12.gvcf.gz
bcftools filter -e 'QUAL<20' --threads 40 -O z geno_AB_snp_an12.gvcf.gz > geno_AB_snp_an12_ql20.gvcf.gz
zcat geno_AB_snp_an12_ql20.gvcf.gz | awk '{if ($5 !~ /*/) print $0}' | bgzip -c > geno_AB_snp_an12_ql20_star.gvcf.gz
bcftools filter -e 'FMT/GT="het"' --set-GTs . --threads 40 -O z geno_AB_snp_an12_ql20_star.gvcf.gz > geno_AB_snp_an12_ql20_star_hetgt.gvcf.gz
bcftools filter -e 'FMT/DP<10' --set-GTs . --threads 40 -O z geno_AB_snp_an12_ql20_star_hetgt.gvcf.gz > geno_AB_snp_an12_ql20_star_hetgt_dp10.gvcf.gz
tabix -p vcf geno_AB_snp_an12_ql20_star_hetgt_dp10.gvcf.gz


#apply variants to the collapsed assembly
for i in $(cat list); do bcftools consensus -f ../../hap0/Eiffel_hap0.fa -M N -a N -I -s $i -p $i -o consensus_hetgt_soft/$i.fasta geno_AB_snp_an12_ql20_star_hetgt_dp10.gvcf.gz; done


#convert consensus genome alignments to non-interleaved format and split by chromosome
cat *.fasta > Ppr.fasta
cat Ppr.fasta | awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' > Ppr_ni.fasta
grep "chr1\b" -A 1 --no-group-separator Ppr_ni.fasta > Ppr_chr1.fasta
sed -i 's/chr1//g' Ppr_chr1.fasta
grep "chr2\b" -A 1 --no-group-separator Ppr_ni.fasta > Ppr_chr2.fasta
sed -i 's/chr2//g' Ppr_chr2.fasta
grep "chr3\b" -A 1 --no-group-separator Ppr_ni.fasta > Ppr_chr3.fasta
sed -i 's/chr3//g' Ppr_chr3.fasta
grep "chr4\b" -A 1 --no-group-separator Ppr_ni.fasta > Ppr_chr4.fasta
sed -i 's/chr4//g' Ppr_chr4.fasta
grep "chr5\b" -A 1 --no-group-separator Ppr_ni.fasta > Ppr_chr5.fasta
sed -i 's/chr5//g' Ppr_chr5.fasta
grep "chr6\b" -A 1 --no-group-separator Ppr_ni.fasta > Ppr_chr6.fasta
sed -i 's/chr6//g' Ppr_chr6.fasta
grep "chr7\b" -A 1 --no-group-separator Ppr_ni.fasta > Ppr_chr7.fasta
sed -i 's/chr7//g' Ppr_chr7.fasta
grep "chr8\b" -A 1 --no-group-separator Ppr_ni.fasta > Ppr_chr8.fasta
sed -i 's/chr8//g' Ppr_chr8.fasta
grep "chr9\b" -A 1 --no-group-separator Ppr_ni.fasta > Ppr_chr9.fasta
sed -i 's/chr9//g' Ppr_chr9.fasta


#split haplotypic regions into 1 Mb bp bins
for i in $(cat list); do msa_split $i.fasta --windows 1000000,0 --out-root $i; done


#calculate trees per 1 Mio bp bin of longest haplotypic blocks
for i in *.fa; do iqtree2 -nt 40 -s $i -m MFP -B 1000 --redo; done



