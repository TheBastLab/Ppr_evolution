# Horizontal gene transfer (HGT) 
## We used the pipeline from https://github.com/reubwn/hgt
## 1.software version
    diamond  v0.9.21.122
## 2.command line      
    export pep=hap0_UTR_longest.pep
    export gff=hap0_add_UTR.gff3
    
    diamond blastp --sensitive --index-chunks 1 -k 500 -e 1e-5 -p 60 -q $pep -d uniref90.dmnd -a Ppr_diamond_blastp_results
    
    cat <(zcat uniref90.fasta.taxlist.gz) <(diamond view -a Ppr_diamond_blastp_results.daa)|perl -lane 'if(@F==2){$tax{$F[0]}=$F[1];}else{if(exists($tax{$F[1]})){print join("\t",@F,$tax{$F[1]});} else {print join("\t",@F,"NA");}}'| gzip > diamond_results.daa.taxid.gz
    
    diamond_to_HGT_candidates.pl -i diamond_results.daa.taxid.gz -f $pep -p taxdump
    
    HGT_candidates_to_fasta.pl -i diamond_results.daa.taxid.gz -c diamond_results.daa.taxid.gz.HGT_candidates.Metazoa.hU30.CHS90.txt -u uniref90.fasta -f $pep -p taxdump
    
    mkdir mafft_alns
    cd outdir
    
    for f in *.fasta; do echo $f; mafft --auto --quiet --thread 8 $f > ../mafft_alns/${f}.mafft; done
    cd ../mafft_alns/
    mkdir processed_files
    
    COUNT=1;
    
    ##iqtree commands
    
    for file in *mafft;
       do echo $COUNT: $file;
       iqtree-omp -s $file -st AA -nt 16 -quiet -bb 1000 -m TESTNEW -msub nuclear
       mv $file processed_files/;
       COUNT=$[COUNT+1];
    done
    
    mkdir iqtree treefiles
    mv *treefile treefiles/
    mv *bionj *gz *contree *iqtree *log *mldist *model *nex iqtree/
    
    get_location_of_HGT_candidates.pl -i diamond_results.daa.taxid.gz.HGT_results.Metazoa.txt -g $gff -f $pep
