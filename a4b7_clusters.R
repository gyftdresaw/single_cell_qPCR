
# look at initial clusters formed from single cell a4b7+ populations
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

## define a4b7+ clusters ##
a4b7_A.cell_ids = c(55,72,58)
a4b7_B1.cell_ids = c(82,79,57,92,71)
a4b7_B2.cell_ids = c(70,68,67)
a4b7_C.cell_ids = c(83,81,59,56,84,60,91,93)
all_a4b7_clusters.cell_ids = c(a4b7_A.cell_ids,a4b7_B1.cell_ids,a4b7_B2.cell_ids,a4b7_C.cell_ids)
all_a4b7_clusters.labels = c(rep('a4b7_A',length(a4b7_A.cell_ids)),
                             rep('a4b7_B1',length(a4b7_B1.cell_ids)),
                             rep('a4b7_B2',length(a4b7_B2.cell_ids)),
                             rep('a4b7_C',length(a4b7_C.cell_ids)))

# plots of these a4b7 clusters
a4b7_clusters = all_fn_rep[,match(all_a4b7_clusters.cell_ids,as.numeric(all_filter.cell_ids))]
colnames(a4b7_clusters) = paste(all_a4b7_clusters.cell_ids,all_a4b7_clusters.labels)

# a4b7 clusters with NAs not replaced
a4b7_clusters_norep = all_fn[,match(all_a4b7_clusters.cell_ids,as.numeric(all_filter.cell_ids))]
colnames(a4b7_clusters_norep) = paste(all_a4b7_clusters.cell_ids,all_a4b7_clusters.labels)

require(gplots)
# recluster to see how well groups are recapitulated
png('cluster_plots/a4b7_reclustered.png',width=1000,height=1000)
heatmap.2(as.matrix(a4b7_clusters),trace='none',density.info='none',scale='none',
          col=redgreen(100),srtCol=90,cexRow=1.0,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# recluster to see how well groups are recapitulated
png('cluster_plots/a4b7_reclustered_scale.png',width=1000,height=1000)
heatmap.2(as.matrix(a4b7_clusters),trace='none',density.info='none',scale='row',
          col=redgreen(100),breaks=seq(-3,3,length.out=101),
          srtCol=90,cexRow=1.0,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# plot given clusters obtained from complete data
png('cluster_plots/a4b7_clusters_scale.png',width=1000,height=1000)
heatmap.2(as.matrix(a4b7_clusters),trace='none',density.info='none',scale='row',Colv=F,
          col=redgreen(100),breaks=seq(-3,3,length.out=101),
          srtCol=90,cexRow=1.0,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

png('cluster_plots/a4b7_clusters.png',width=1000,height=1000)
heatmap.2(as.matrix(a4b7_clusters),trace='none',density.info='none',scale='none',Colv=F,
          col=redgreen(100),srtCol=90,cexRow=1.0,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# correlation heatmaps
require(RColorBrewer)
bluescale = colorRampPalette(brewer.pal(9,"Blues"))(100)
quakescale = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")

# across cells
png('cluster_plots/a4b7_clusters_cor_heatmap.png',width=1000,height=1000)
heatmap.2(cor(as.matrix(a4b7_clusters)),Colv='Rowv',trace='none',density.info='none',col=bluescale,
          srtCol=90,cexRow=1.5,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# across genes
a4b7_gene_sds = apply(a4b7_clusters,1,sd)
a4b7_clusters_var = a4b7_clusters[a4b7_gene_sds > 0,]

png('cluster_plots/a4b7_clusters_gg_cor_heatmap.png',width=1200,height=1200)
heatmap.2(cor(as.matrix(t(a4b7_clusters_var))),Colv='Rowv',trace='none',density.info='none',col=bluescale,
          srtCol=90,cexRow=1.5,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# sorting
# tf_cor = cor(as.matrix(t(a4b7_clusters_var)))
# sort(tf_cor[rownames(tf_cor)=='Zbtb16',])

# PCA of different clusters
a4b7_clusters_pca = prcomp(t(a4b7_clusters_var),scale=T)

a4b7_clusters_pca.df = data.frame(a4b7_clusters_pca$x[,1:10])
a4b7_clusters_pca.df$Type = all_a4b7_clusters.labels

var_explained = (a4b7_clusters_pca$sdev)^2 / sum(a4b7_clusters_pca$sdev^2)

g = ggplot(a4b7_clusters_pca.df, aes(x=PC1,y=PC2,color=Type)) + geom_point(size=4)
g = g + xlab('PC1 (16.46%)') + ylab('PC2 (10.69%)')
g = g + theme(axis.text=element_text(size=16))
g = g + theme(axis.title=element_text(size=16,face='bold'))
g

ggsave('cluster_plots/a4b7_clusters_pca12.png')

g = ggplot(a4b7_clusters_pca.df, aes(x=PC2,y=PC3,color=Type)) + geom_point(size=4)
g = g + xlab('PC2 (10.69%)') + ylab('PC3 (9.03%)')
g = g + theme(axis.text=element_text(size=16))
g = g + theme(axis.title=element_text(size=16,face='bold'))
g

ggsave('cluster_plots/a4b7_clusters_pca23.png')

# sort(abs(a4b7_clusters_pca$rotation[,1]))

## want to find genes distinguishing different clusters ##

# linear fit over groups on all genes to find most distinguishing genes
a4b7_linfits = vector('list',dim(a4b7_clusters)[1])
a4b7_anova = vector('list',dim(a4b7_clusters)[1])
for(i in 1:dim(a4b7_clusters)[1]){
  gene_index = i
  anova_test.df = data.frame(level=as.numeric(a4b7_clusters[gene_index,]),cluster=all_a4b7_clusters.labels)
  
  if(sd(anova_test.df$level) > 0){
    linfit = lm(level~cluster,anova_test.df)
    a4b7_linfits[[i]] = linfit
    a4b7_anova[[i]] = anova(linfit)
  }
  # for printing specific results
  # summary(linfit)
  # anova(linfit)
  # rownames(ilcp_clusters)[gene_index]
}

# extract F test p values
a4b7_anova_fvals = sapply(a4b7_anova,
                          function(x) if(is.null(x)) {return(NA)} else{return(x[[5]][1])})
a4b7_anova_ind = a4b7_anova_fvals < 0.10
a4b7_anova_ind[is.na(a4b7_anova_ind)] = F
# add some bonus genes
a4b7_anova_ind[match('Zbtb16',rownames(a4b7_clusters))] = T
a4b7_anova_ind[match('Bcl11a',rownames(a4b7_clusters))] = T
a4b7_anova_genes = rownames(a4b7_clusters)[a4b7_anova_ind]
a4b7_anova_genes

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
a4b7_splot = function(gene,save=F){
  gene_index = match(gene,rownames(a4b7_clusters))
  anova_test.df = data.frame(level=as.numeric(a4b7_clusters[gene_index,]),cluster=all_a4b7_clusters.labels)
  
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
    ggsave(paste('cluster_plots/',gene,'_levels_a4b7.png',sep=''))
  }
}

# only 17 out of 93 genes significant at 0.10
# length(a4b7_anova_genes)

# heatmap these specific genes
png('cluster_plots/a4b7_clusters_anova.png',width=1000,height=1000)
heatmap.2(as.matrix(a4b7_clusters[a4b7_anova_ind,]),trace='none',density.info='none',scale='none',Colv=F,
          col=redgreen(100),
          srtCol=90,cexRow=1.7,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# row scale
png('cluster_plots/a4b7_clusters_anova_scale.png',width=1000,height=1000)
heatmap.2(as.matrix(a4b7_clusters[a4b7_anova_ind,]),trace='none',density.info='none',scale='row',Colv=F,
          col=redgreen(100),breaks=seq(-3,3,length.out=101),
          srtCol=90,cexRow=1.7,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# heatmap specific genes with reclustering
png('cluster_plots/a4b7_clusters_anova_recluster.png',width=1000,height=1000)
heatmap.2(as.matrix(a4b7_clusters[a4b7_anova_ind,]),trace='none',density.info='none',scale='none',Colv=T,
          col=redgreen(100),
          srtCol=90,cexRow=1.7,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# row scale recluster
png('cluster_plots/a4b7_clusters_anova_scale_recluster.png',width=1000,height=1000)
heatmap.2(as.matrix(a4b7_clusters[a4b7_anova_ind,]),trace='none',density.info='none',scale='row',Colv=T,
          col=redgreen(100),breaks=seq(-3,3,length.out=101),
          srtCol=90,cexRow=1.7,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# average across cells in groups
agg_result = aggregate(t(a4b7_clusters),by=list(ID=all_a4b7_clusters.labels),mean)
agg_a4b7_clusters = t(agg_result[,2:dim(agg_result)[2]])
colnames(agg_a4b7_clusters) = agg_result$ID

# heatmap average results
png('cluster_plots/a4b7_clusters_ave_anova.png',width=1000,height=1000)
heatmap.2(as.matrix(agg_a4b7_clusters[a4b7_anova_ind,]),trace='none',density.info='none',scale='none',Colv=T,
          col=redgreen(100),breaks=seq(-30,0,length.out=101),
          srtCol=90,cexRow=1.7,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# row scaled
png('cluster_plots/a4b7_clusters_ave_anova_scale.png',width=1000,height=1000)
heatmap.2(as.matrix(agg_a4b7_clusters[a4b7_anova_ind,]),trace='none',density.info='none',scale='row',Colv=T,
          col=redgreen(100),srtCol=90,cexRow=1.7,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

