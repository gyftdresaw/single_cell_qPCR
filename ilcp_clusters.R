
# look at initial clusters formed from single cell ILCP
# hierarchical clusters based on all samples isolated

setwd('~/Documents/Dinner/computer transfer/Bendelac/single-cell-qPCR/')

# load combined threshold adjusted data
load('data/parsed/adj_combined.RData')

# need to normalize data by housekeeping genes
# should be all set now, column names are the same as chip1
hkg_rows = all_chips.df$Symbol %in% c('Hprt','Gapdh','Actb')
hkg_pass = apply(all_chips.df[hkg_rows,-1], 2, function(x) !any(is.na(x)))

# getting rid of first symbol column
all_filter = all_chips.df[,c(F,hkg_pass)]
rownames(all_filter) = all_chips.df[,1]
all_filter.types = all_chips.cell_types[c(F,hkg_pass)]
all_filter.cell_ids = all_chips.cell_ids[c(F,hkg_pass)]
colnames(all_filter) = paste(all_filter.cell_ids,all_filter.types)

# before normalization, replaced
all_filter_rep = all_filter
all_filter_rep[is.na(all_filter_rep)] = 40
all_filter_rep = -all_filter_rep

# normalize by housekeeping genes
hkg_ave = apply(all_filter[hkg_rows,], 2, mean)
# normalized results are -deltaCt
all_fn = -1 * sweep(all_filter,MARGIN=2,hkg_ave,FUN="-")

all_fn_rep = all_fn
all_fn_rep[is.na(all_fn_rep)] = -30

### define ILCP clusters ###
ilcp_A.cell_ids = c(54,49,51,62,75,77,65,86)
ilcp_B.cell_ids = c(52,73,63,64)
ilcp_C.cell_ids = c(50,53,85,76,89,88)
ilcp_D.cell_ids = c(66,74,87)
all_ilcp_clusters.cell_ids = c(ilcp_A.cell_ids,ilcp_B.cell_ids,ilcp_C.cell_ids,ilcp_D.cell_ids)
all_ilcp_clusters.labels = c(rep('ILCP_A',length(ilcp_A.cell_ids)),
                             rep('ILCP_B',length(ilcp_B.cell_ids)),
                             rep('ILCP_C',length(ilcp_C.cell_ids)),
                             rep('ILCP_D',length(ilcp_D.cell_ids)))

# plots of these ilcp clusters
ilcp_clusters = all_fn_rep[,match(all_ilcp_clusters.cell_ids,as.numeric(all_filter.cell_ids))]
colnames(ilcp_clusters) = paste(all_ilcp_clusters.cell_ids,all_ilcp_clusters.labels)

# ilcp clusters with NAs not replaced
ilcp_clusters_norep = all_fn[,match(all_ilcp_clusters.cell_ids,as.numeric(all_filter.cell_ids))]
colnames(ilcp_clusters_norep) = paste(all_ilcp_clusters.cell_ids,all_ilcp_clusters.labels)

require(gplots)
# recluster to see how well groups are recapitulated
png('cluster_plots/ilcp_reclustered.png',width=1000,height=1000)
heatmap.2(as.matrix(ilcp_clusters),trace='none',density.info='none',scale='none',
          col=redgreen(100),srtCol=90,cexRow=1.0,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# recluster to see how well groups are recapitulated
png('cluster_plots/ilcp_reclustered_scale.png',width=1000,height=1000)
heatmap.2(as.matrix(ilcp_clusters),trace='none',density.info='none',scale='row',
          col=redgreen(100),breaks=seq(-3,3,length.out=101),
          srtCol=90,cexRow=1.0,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# plot given clusters obtained from complete data
png('cluster_plots/ilcp_clusters_scale.png',width=1000,height=1000)
heatmap.2(as.matrix(ilcp_clusters),trace='none',density.info='none',scale='row',Colv=F,
          col=redgreen(100),breaks=seq(-3,3,length.out=101),
          srtCol=90,cexRow=1.0,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

png('cluster_plots/ilcp_clusters.png',width=1000,height=1000)
heatmap.2(as.matrix(ilcp_clusters),trace='none',density.info='none',scale='none',Colv=F,
          col=redgreen(100),srtCol=90,cexRow=1.0,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# correlation heatmaps
require(RColorBrewer)
bluescale = colorRampPalette(brewer.pal(9,"Blues"))(100)
quakescale = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")

# across cells
png('cluster_plots/ilcp_clusters_cor_heatmap.png',width=1000,height=1000)
heatmap.2(cor(as.matrix(ilcp_clusters)),Colv='Rowv',trace='none',density.info='none',col=bluescale,
          srtCol=90,cexRow=1.5,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# across genes
ilcp_gene_sds = apply(ilcp_clusters,1,sd)
ilcp_clusters_var = ilcp_clusters[ilcp_gene_sds > 0,]

png('cluster_plots/ilcp_clusters_gg_cor_heatmap.png',width=1200,height=1200)
heatmap.2(cor(as.matrix(t(ilcp_clusters_var))),Colv='Rowv',trace='none',density.info='none',col=bluescale,
          srtCol=90,cexRow=1.5,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# sorting
# tf_cor = cor(as.matrix(t(ilcp_clusters_var)))
# sort(tf_cor[rownames(tf_cor)=='Zbtb16',])

# PCA of different clusters
ilcp_gene_sds = apply(ilcp_clusters,1,sd)
ilcp_clusters_var = ilcp_clusters[ilcp_gene_sds > 0,]
ilcp_clusters_pca = prcomp(t(ilcp_clusters_var),scale=T)

ilcp_clusters_pca.df = data.frame(ilcp_clusters_pca$x[,1:10])
ilcp_clusters_pca.df$Type = all_ilcp_clusters.labels

var_explained = (ilcp_clusters_pca$sdev)^2 / sum(ilcp_clusters_pca$sdev^2)

g = ggplot(ilcp_clusters_pca.df, aes(x=PC1,y=PC2,color=Type)) + geom_point(size=4)
g = g + xlab('PC1 (17.81%)') + ylab('PC2 (10.84%)')
g = g + theme(axis.text=element_text(size=16))
g = g + theme(axis.title=element_text(size=16,face='bold'))
g

ggsave('cluster_plots/ilcp_clusters_pca12.png')

g = ggplot(ilcp_clusters_pca.df, aes(x=PC1,y=PC3,color=Type)) + geom_point(size=4)
g = g + xlab('PC1 (17.81%)') + ylab('PC3 (9.96%)')
g = g + theme(axis.text=element_text(size=16))
g = g + theme(axis.title=element_text(size=16,face='bold'))
g

ggsave('cluster_plots/ilcp_clusters_pca13.png')

# sort(abs(ilcp_clusters_pca$rotation[,1]))

## want to find genes distinguishing different clusters ##

# linear fit over groups on all genes to find most distinguishing genes
ilcp_linfits = vector('list',dim(ilcp_clusters)[1])
ilcp_anova = vector('list',dim(ilcp_clusters)[1])
for(i in 1:dim(ilcp_clusters)[1]){
  gene_index = i
  anova_test.df = data.frame(level=as.numeric(ilcp_clusters[gene_index,]),cluster=all_ilcp_clusters.labels)
  
  if (sd(anova_test.df$level) > 0){
    linfit = lm(level~cluster,anova_test.df)
    ilcp_linfits[[i]] = linfit
    ilcp_anova[[i]] = anova(linfit)
  }
  # for printing specific results
  # summary(linfit)
  # anova(linfit)
  # rownames(ilcp_clusters)[gene_index]
}

# extract F test p values
ilcp_anova_fvals = sapply(ilcp_anova,
                          function(x) if(is.null(x)) {return(NA)} else {return(x[[5]][1])})
ilcp_anova_ind = ilcp_anova_fvals < 0.05
ilcp_anova_ind[is.na(ilcp_anova_ind)] = F
ilcp_anova_genes = rownames(ilcp_clusters)[ilcp_anova_ind]

# summarySE helper function
#############################
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#############
require(ggplot2)
# make scatterplot of individual gene
ilcp_splot = function(gene,save=F){
  gene_index = match(gene,rownames(ilcp_clusters))
  anova_test.df = data.frame(level=as.numeric(ilcp_clusters[gene_index,]),cluster=all_ilcp_clusters.labels)
  
  # use helper function found below
  anova_sum = summarySE(anova_test.df,measurevar="level",groupvars=c("cluster"))
  
  g = ggplot(anova_test.df,aes(x=cluster,y=level))
  g = g + geom_point(size=4,position=position_jitter(width=0.3))
  g = g + ggtitle(rownames(ilcp_clusters)[gene_index])
  g = g + theme(axis.title=element_text(size=16,face='bold'))
  g = g + theme(axis.text=element_text(size=16))
  g = g + theme(title=element_text(size=16))
  print(g)
  
  if (save){
    ggsave(paste('cluster_plots/',gene,'_levels_ilcp.png',sep=''))
  }
}

# 31 genes out of 93 with p value < 0.05
# sum(anova_fvals < 0.05)

# heatmap these specific genes
png('cluster_plots/ilcp_clusters_anova.png',width=1000,height=1000)
heatmap.2(as.matrix(ilcp_clusters[ilcp_anova_ind,]),trace='none',density.info='none',scale='none',Colv=F,
          col=redgreen(100),
          srtCol=90,cexRow=1.7,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# row scale
png('cluster_plots/ilcp_clusters_anova_scale.png',width=1000,height=1000)
heatmap.2(as.matrix(ilcp_clusters[ilcp_anova_ind,]),trace='none',density.info='none',scale='row',Colv=F,
          col=redgreen(100),breaks=seq(-3,3,length.out=101),
          srtCol=90,cexRow=1.7,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# heatmap specific genes with reclustering
png('cluster_plots/ilcp_clusters_anova_recluster.png',width=1000,height=1000)
heatmap.2(as.matrix(ilcp_clusters[ilcp_anova_ind,]),trace='none',density.info='none',scale='none',Colv=T,
          col=redgreen(100),
          srtCol=90,cexRow=1.7,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# row scale recluster
png('cluster_plots/ilcp_clusters_anova_scale_recluster.png',width=1000,height=1000)
heatmap.2(as.matrix(ilcp_clusters[ilcp_anova_ind,]),trace='none',density.info='none',scale='row',Colv=T,
          col=redgreen(100),breaks=seq(-3,3,length.out=101),
          srtCol=90,cexRow=1.7,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# average across cells in groups
agg_result = aggregate(t(ilcp_clusters),by=list(ID=all_ilcp_clusters.labels),mean)
agg_ilcp_clusters = t(agg_result[,2:dim(agg_result)[2]])
colnames(agg_ilcp_clusters) = agg_result$ID

# heatmap average results
png('cluster_plots/ilcp_clusters_ave_anova.png',width=1000,height=1000)
heatmap.2(as.matrix(agg_ilcp_clusters[ilcp_anova_ind,]),trace='none',density.info='none',scale='none',Colv=T,
          col=redgreen(100),breaks=seq(-30,0,length.out=101),
          srtCol=90,cexRow=1.7,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# row scaled
png('cluster_plots/ilcp_clusters_ave_anova_scale.png',width=1000,height=1000)
heatmap.2(as.matrix(agg_ilcp_clusters[ilcp_anova_ind,]),trace='none',density.info='none',scale='row',Colv=T,
          col=redgreen(100),srtCol=90,cexRow=1.7,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

################################
# fishers/chisquare testing of binary expression 
#
# use data without NA replaced

# linear fit over groups on all genes to find most distinguishing genes
ilcp_fishtests = vector('list',dim(ilcp_clusters_norep)[1])
ilcp_chitests = vector('list',dim(ilcp_clusters_norep)[1])
ilcp_norep_linfits = vector('list',dim(ilcp_clusters_norep)[1])
ilcp_norep_anova = vector('list',dim(ilcp_clusters_norep)[1])
for(i in 1:dim(ilcp_clusters_norep)[1]){
  gene_index = i
  bin_test.df = data.frame(level=as.numeric(ilcp_clusters_norep[gene_index,]),cluster=all_ilcp_clusters.labels)
  
  bin_table = table(factor(!is.na(bin_test.df$level),levels=c(FALSE,TRUE)),bin_test.df$cluster)
  
  # only try fisher test if there are mixed proportions
  if (sum(rowSums(bin_table) == 0) == 0){
    ilcp_fishtests[[i]] = fisher.test(bin_table)
    ilcp_chitests[[i]] = chisq.test(bin_table,simulate.p.value=T)
  }
  
  # only try to fit linear model if more than one cluster is represented
  if (length(unique(bin_test.df[!is.na(bin_test.df$level),]$cluster)) > 1) {
    ilcp_norep_linfits[[i]] = lm(level~cluster,bin_test.df)
    ilcp_norep_anova[[i]] = anova(ilcp_norep_linfits[[i]])
  }
}

# extract p values
fisher_pvals = sapply(ilcp_fishtests,
                      function(x) if (is.null(x)) {return(NA)} else{return(x[[1]])})
chisq_pvals = sapply(ilcp_chitests,
                     function(x) if (is.null(x)) {return(NA)} else{return(x[[3]])})
# need to check for nulls (incomplete categories for anova)
norep_anova_fvals = sapply(ilcp_norep_anova,
                           function(x) if (is.null(x)){return(NA)}else{return(x[[5]][1])})

# these are pretty much the same, good correlation
# going to just go with fisher test
# fisher test also correlates well with simple anova from above
fisher_genes = rownames(ilcp_clusters_norep)[fisher_pvals < 0.05]
# remove na vals
fisher_genes = fisher_genes[!is.na(fisher_genes)]
norep_anova_genes = rownames(ilcp_clusters_norep)[norep_anova_fvals < 0.05]
# remove na vals
norep_anova_genes = norep_anova_genes[!is.na(norep_anova_genes)]

# concatenate all binary test based genes
# 35 genes in total
bin_genes = union(fisher_genes,norep_anova_genes)

# previously 31 genes were found just by anova
# 30 are shared between previous anova and binary test strategy
# 
# Tbp, Runx1, Tet2, Cxcr6, Il-7r uniquely found by binary test
# Notch1 uniquely found by previous anova

# prestring to label genes (rownames) by significance in fisher or norep anova
fish_pre = rep('',dim(ilcp_clusters_norep)[1])
fish_pre[fisher_pvals < 0.05] = '#'
norep_anova_pre = rep('',dim(ilcp_clusters_norep)[1])
norep_anova_pre[norep_anova_fvals < 0.05] = '*'

sig_gene_labels = paste(fish_pre,norep_anova_pre,rownames(ilcp_clusters_norep),sep='')

# heatmap genes found by binary test with significance labels
# first heatmap over all samples

# logical indicator of genes of interest
bin_ind = norep_anova_fvals < 0.05
bin_ind[is.na(bin_ind)] = FALSE
fish_ind = fisher_pvals < 0.05
fish_ind[is.na(fish_ind)] = FALSE
bin_ind = bin_ind | fish_ind

# set new rownames
bin_ilcp_clusters = ilcp_clusters
rownames(bin_ilcp_clusters) = sig_gene_labels

# heatmap these specific genes
png('cluster_plots/ilcp_clusters_bin.png',width=1000,height=1000)
heatmap.2(as.matrix(bin_ilcp_clusters[bin_ind,]),trace='none',density.info='none',scale='none',Colv=F,
          col=redgreen(100),
          srtCol=90,cexRow=1.7,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# row scale
png('cluster_plots/ilcp_clusters_bin_scale.png',width=1000,height=1000)
heatmap.2(as.matrix(bin_ilcp_clusters[bin_ind,]),trace='none',density.info='none',scale='row',Colv=F,
          col=redgreen(100),breaks=seq(-3,3,length.out=101),
          srtCol=90,cexRow=1.7,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# heatmap specific binary genes wtih reclustering
png('cluster_plots/ilcp_clusters_bin_recluster.png',width=1000,height=1000)
heatmap.2(as.matrix(bin_ilcp_clusters[bin_ind,]),trace='none',density.info='none',scale='none',Colv=T,
          col=redgreen(100),
          srtCol=90,cexRow=1.7,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# row scale reclustering
png('cluster_plots/ilcp_clusters_bin_scale_recluster.png',width=1000,height=1000)
heatmap.2(as.matrix(bin_ilcp_clusters[bin_ind,]),trace='none',density.info='none',scale='row',Colv=T,
          col=redgreen(100),breaks=seq(-3,3,length.out=101),
          srtCol=90,cexRow=1.7,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# averages, set rownames
bin_agg_ilcp_clusters = agg_ilcp_clusters
rownames(bin_agg_ilcp_clusters) = sig_gene_labels

# heatmap average results
png('cluster_plots/ilcp_clusters_ave_bin.png',width=1000,height=1000)
heatmap.2(as.matrix(bin_agg_ilcp_clusters[bin_ind,]),trace='none',density.info='none',scale='none',Colv=T,
          col=redgreen(100),breaks=seq(-30,0,length.out=101),
          srtCol=90,cexRow=1.7,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# row scaled
png('cluster_plots/ilcp_clusters_ave_bin_scale.png',width=1000,height=1000)
heatmap.2(as.matrix(bin_agg_ilcp_clusters[bin_ind,]),trace='none',density.info='none',scale='row',Colv=T,
          col=redgreen(100),srtCol=90,cexRow=1.7,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()


