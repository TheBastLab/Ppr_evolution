# *Platynothrus peltifer* genome analysis
*Platynothrus peltifer* genome with haplotype-specific analyses 

## Table of contents
* [Initial analyses on HiFi reads](https://github.com/TheBastLab/Ppr_evolution/blob/main/initial_hifi_analyses.md)
* [Collapsed chromosome-level assembly pipeline for German Ppr](https://github.com/TheBastLab/Ppr_evolution/blob/main/assembly_hap0_German_Ppr.md)
* [Phased assembly pipeline](https://github.com/TheBastLab/Ppr_evolution/blob/main/phased_assembly.md)
* [Assembly evaluation](https://github.com/TheBastLab/Ppr_evolution/blob/main/assembly_evaluation.md)
* [Variant calling](https://github.com/TheBastLab/Ppr_evolution/blob/main/variant_calling.md)
* [Genome annotation](https://github.com/TheBastLab/Ppr_evolution/blob/main/Genome_annotation.md)
* [Population statistics](https://github.com/TheBastLab/Ppr_evolution/blob/main/population_statistics.md)
* [Differential expression of alleles](https://github.com/TheBastLab/Ppr_evolution/blob/main/Differential_expression_of_alleles.md)
* [Horizontal gene transfer](https://github.com/TheBastLab/Ppr_evolution/blob/main/Horizontal_gene_transfer.md)
* [Wolbachia](https://github.com/TheBastLab/Ppr_evolution/blob/main/Wolbachia.md)

Steps followed to perform MDS analyses and explore evolutionary dynamics of the two diverging haplotypes (alleles) in an asexually reproducing mite:

1. **MDS Analyses**: For detailed steps on performing MDS analyses, please check [MDS](https://github.com/TheBastLab/Ppr_evolution/blob/main/MDSing.md).
2. **Diversity Analysis**: For information on overall as well as nucleotide diversity at synonymous/nonsynonymous sites, please check [here](https://github.com/TheBastLab/Ppr_evolution/blob/main/calculateDiversity.md).
3. **Omega Calculation**: To calculate omega using HyPhy, please check the [evolutionary rates pipeline](https://github.com/TheBastLab/Ppr_evolution/blob/main/calculate_omega.md).
4. **SFS Profiles**: Please use the following [script](https://github.com/merrbii/Ppr_hap_divergence/blob/main/src/get_sfs_profiles.sh) that processes a VCF file to extract, count, and sort unique genotypes. It then formats the results and saves them to an output file that can be used to produce bar plots in R, for instance.
