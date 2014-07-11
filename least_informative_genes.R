
# identify least informative genes for clustering 
# in order to modify qPCR assay to be ILC specific

# based on data/vars in initial.R

# 10 cell and single cell counts

ss_na_counts = apply(!is.na(all_fn[,bulkss_types=='Single Cell']),1,sum)
bulk_na_counts = apply(!is.na(all_fn[,bulkss_types=='10 Cell']),1,sum)

# lowest PCA loadings
# based on only single cell PCA

N = 3 # number of top PC's
sort(apply(abs(ss_single_pca$rotation[,1:N]),1,max))
