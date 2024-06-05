## Calculating MDS Plots Using PLINK

This pipeline describes the steps to perform MDS using PLINK, starting from a VCF file.

### Step 1: Convert VCF to PLINK format

```bash
plink2 --allow-extra-chr --double-id --make-bed --max-alleles 2 --min-alleles 2 --out pca/hapAB.allPops.miss0.8.var.ann.fixedGT.biall --set-missing-var-ids @:# --vcf hapAB.allPops.miss0.8.var.ann.fixedGT.vcf.gz
```

This command converts the VCF file to PLINK binary format, ensuring only biallelic sites are included.

### Step 2: Perform LD pruning

```bash
plink2 --allow-extra-chr --bfile hapAB.allPops.miss0.8.var.ann.fixedGT.biall --indep-pairwise 1 kb 1 0.2 --out LDpruned/hapAB.allPops.miss0.8.var.ann.fixedGT.biall.ldkb1r0.2
```

This step performs linkage disequilibrium (LD) pruning with a window size of 1 kb, step size of 1, and an r-squared threshold of 0.2.

### Step 3: Create pruned PLINK files

```bash
plink2 --allow-extra-chr --bfile hapAB.allPops.miss0.8.var.ann.fixedGT.biall --exclude LDpruned/hapAB.allPops.miss0.8.var.ann.fixedGT.biall.ldkb1r0.2.prune.out --make-bed --out LDpruned/hapAB.allPops.miss0.8.var.ann.fixedGT.biall.ldkb1r0.2.pruned
```

This command creates new PLINK binary files excluding the SNPs identified during LD pruning.

### Step 4: Perform MDS clustering
```bash
plink --allow-extra-chr --bfile hapAB.allPops.miss0.8.var.ann.fixedGT.biall.ldkb1r0.2.pruned --cluster --mds-plot 5 eigvals --out hapAB.allPops.miss0.8.var.ann.fixedGT.biall.ldkb1r0.2.pruned.mds
```

This step performs MDS clustering and generates a plot with the top 5 eigenvalues.
