# Wolbachia
## 1.Software version
diamond  v0.9.21.122
NGenomeSyn  v.141

## 2.Commandlines

    diamond blastp --sensitive --index-chunks 1 -k 500 -e 1e-5 -p 60 -q GCF_016584425.1_ASM1658442v1_protein.faa -d Ppr.pep.dmnd -a ASM1658442v1_Ppr_hap0
    diamond view -a ASM1658442v1_Ppr_hap0.daa -o ASM1658442v1_Ppr_hap0.blast

    NGenomeSyn -InConf ASM1658442v1_Ppr_hap0.conf -OutPut ASM1658442v1_Ppr_hap0

