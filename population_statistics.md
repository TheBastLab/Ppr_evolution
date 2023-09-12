## Population Statistics

For steps 1~4, each individual is processed separately by an independent run of the Perl script; shown here is DE_T501 as an example.
For step 5, each pair of individual (within population) is processed separately by an independent run of the Perl script; shown here is DE_T501 to DE_T502 as an example.

### Step 1: Obtain coverage distribution 
```
perl complete_5pop_popstats_1_coverage_distribution.pl DE_T501
```
### Step 1.5: Calculate coverage threshod (in R)
```
d0 = read.table(paste("DE_T501","_coverage_distribution.txt",sep=""),header=TRUE)
Cov_median = median(rep(d0[,1],d0[,2]))
Cov_thres = ceiling(Cov_median*0.75)
```
### Step 2: Produce coverage+repeat genomic masking information
```
perl complete_5pop_popstats_2_coverage_labelling.pl DE_T501 73
```
### Step 3: Calculate individual heterozygosity by 1MB blocks
```
perl complete_5pop_popstats_3_individual_heterozygosity.pl DE_T501
```
### Step 4: Calculate distribution of homozygote-stretch lengths
```
perl complete_5pop_popstats_4_homozygote_stretch.pl DE_T501
```
### Step 5: Calculate shared heterozygosity between individual pairs
```
perl complete_5pop_popstats_5_pairwise_shared_het.pl DE_T501 DE_T502
```
