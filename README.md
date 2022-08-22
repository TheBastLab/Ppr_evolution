# Ppr_evolution
Ppr genome with haplotpye specific analyses and populations



This is the github for the Ppr genome project.


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
