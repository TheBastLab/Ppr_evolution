## Initial analyses on HiFi reads

### *k*-mer analysis

[KAT](https://github.com/TGAC/KAT) version 2.4.2
```sh
kat hist -o kat_hist hifi_reads.fastq.gz
kat gcp -o kat_gcp hifi_reads.fastq.gz
```

<img src="./fig/kat_hist.svg" width=400> <img src="./fig/kat_gcp.svg" width=400>

[Smudgeplot](https://github.com/KamilSJaron/smudgeplot) version 0.2.5
[KMC](https://github.com/tbenavi1/KMC) version 
```sh
mkdir tmp_smudge
ls hifi_reads.fastq.gz > FILES
kmc -k27 -ci1 -cs10000 @FILES kmcdb tmp_smudge
kmc_tools transform kmcdb histogram kmcdb_k27.hist -cx10000

L=$(smudgeplot.py cutoff kmcdb_k27.hist L)
U=$(smudgeplot.py cutoff kmcdb_k27.hist U)
echo $L $U

kmc_tools transform kmcdb -ci"$L" -cx"$U" dump -s kmcdb_L"$L"_U"$U".dump
smudgeplot.py hetkmers -o kmcdb_L"$L"_U"$U" < kmcdb_L"$L"_U"$U".dump

smudgeplot.py plot kmcdb_L"$L"_U"$U"_coverages.tsv
```
