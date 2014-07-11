
# initial looks at the data
require(ggplot2)

setwd('~/Documents/Dinner/computer transfer/Bendelac/single-cell-qPCR/')

combined_file = 'data/parsed/combined_clean.txt'
# read first two header lines separately 
qPCR_load = function(data_file){
  cell_ids = scan(data_file,nlines=1,what=character(),sep='\t')
  cell_types = scan(data_file,skip=1,nlines=1,what=character(),sep='\t')
  # now read data
  qPCR.df = read.table(data_file,skip=2,sep='\t')
  return(list(qPCR.df,cell_ids,cell_types))
}

qPCR.lst = qPCR_load(combined_file)
qPCR.df = qPCR.lst[[1]]
qPCR.cell_ids = qPCR.lst[[2]]
qPCR.cell_types = qPCR.lst[[3]]

# load up all three chips individually
chip1.lst = qPCR_load('data/parsed/chip1_clean.txt')
chip1.df = chip1.lst[[1]]
chip1.cell_ids = chip1.lst[[2]]
chip1.cell_types = chip1.lst[[3]]
chip1.thresh = 0.096603277 # fluorescence threshold from experiment

chip2.lst = qPCR_load('data/parsed/chip2_clean.txt')
chip2.df = chip2.lst[[1]]
chip2.cell_ids = chip2.lst[[2]]
chip2.cell_types = chip2.lst[[3]]
chip2.thresh = 0.032601407 # fluorescence threshold from experiment

chip3.lst = qPCR_load('data/parsed/chip3_clean.txt')
chip3.df = chip3.lst[[1]]
chip3.cell_ids = chip3.lst[[2]]
chip3.cell_types = chip3.lst[[3]]
chip3.thresh = 0.036386395 # fluorescence threshold from experiment

# plot histograms of all non-NA values
chip1.vals = na.omit(c(as.matrix(chip1.df[,-1])))
chip2.vals = na.omit(c(as.matrix(chip2.df[,-1])))
chip3.vals = na.omit(c(as.matrix(chip3.df[,-1])))

all_vals.df = data.frame(vals=c(chip1.vals,chip2.vals,chip3.vals),
                         chips=c(rep('chip1',length(chip1.vals)),rep('chip2',length(chip2.vals)),
                                 rep('chip3',length(chip3.vals))))

ggplot(all_vals.df) + 
  geom_histogram(data=subset(all_vals.df,chips=='chip1'),aes(x=vals,y=..density..),fill='red',alpha=0.2,position='identity') + 
  geom_histogram(data=subset(all_vals.df,chips=='chip2'),aes(x=vals,y=..density..),fill='blue',alpha=0.2,position='identity') + 
  geom_histogram(data=subset(all_vals.df,chips=='chip3'),aes(x=vals,y=..density..),fill='green',alpha=0.2,position='identity')

ggsave('plots/chip_hist.png')

## adjust by corresponding Ct thresh
all_vals_adj.df = data.frame(vals=c(chip1.vals,chip2.vals+log2(chip1.thresh/chip2.thresh),
                                    chip3.vals+log2(chip1.thresh/chip3.thresh)),
                            chips=c(rep('chip1',length(chip1.vals)),rep('chip2',length(chip2.vals)),
                                    rep('chip3',length(chip3.vals))))

ggplot(all_vals_adj.df) + 
  geom_histogram(data=subset(all_vals_adj.df,chips=='chip1'),aes(x=vals,y=..density..),fill='red',alpha=0.2,position='identity') + 
  geom_histogram(data=subset(all_vals_adj.df,chips=='chip2'),aes(x=vals,y=..density..),fill='blue',alpha=0.2,position='identity') + 
  geom_histogram(data=subset(all_vals_adj.df,chips=='chip3'),aes(x=vals,y=..density..),fill='green',alpha=0.2,position='identity')

ggsave('plots/chip_hist_adj.png')

## going to proceed with the adjustment
chip2.df[,-1] = chip2.df[,-1] + log2(chip1.thresh/chip2.thresh)
chip3.df[,-1] = chip3.df[,-1] + log2(chip1.thresh/chip3.thresh)

# now combine all chips into one data matrix
# colnames given by cell ids
colnames(chip1.df) = c('Symbol',chip1.cell_ids[-1])
colnames(chip2.df) = c('Symbol',chip2.cell_ids[-1])
colnames(chip3.df) = c('Symbol',chip3.cell_ids[-1])

# merge chips 2 and 3 by shared Symbol
chip23.df = merge(chip2.df,chip3.df,by='Symbol')

# now we want to merge with chip1 by column names
all_chips.df = rbind(chip1.df,chip23.df)
all_chips.cell_types = chip1.cell_types
all_chips.cell_ids = chip1.cell_ids
# save(all_chips.df,all_chips.cell_types,all_chips.cell_ids,file='data/parsed/adj_combined.RData')

# should be all set now, column names are the same as chip1
hkg_rows = all_chips.df$Symbol %in% c('Hprt','Gapdh','Actb')
hkg_pass = apply(all_chips.df[hkg_rows,-1], 2, function(x) !any(is.na(x)))

# getting rid of first symbol column
all_filter = all_chips.df[,c(F,hkg_pass)]
rownames(all_filter) = all_chips.df[,1]
all_filter.types = chip1.cell_types[c(F,hkg_pass)]
all_filter.cell_ids = chip1.cell_ids[c(F,hkg_pass)]
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

gene_sds = apply(all_fn_rep,1,sd)
gene_means = apply(all_fn_rep,1,mean)

hkg = c('Hprt','Gapdh','Actb')
cv.df = data.frame(sds = gene_sds,means=gene_means,symbol=rownames(all_fn_rep))
g = ggplot(data=cv.df,mapping=aes(means,sds))
g = g + geom_point(colour='blue',size=4)
g = g + geom_point(mapping=aes(means,sds),data=cv.df[cv.df$symbol %in% hkg,],colour='red',size=4)
g = g + geom_text(mapping=aes(means-0.1,sds+0.2,label=symbol),data=cv.df[cv.df$symbol %in% hkg,],hjust=1,vjust=0,colour='red',size=6,angle=-60)
g = g + xlab('Mean') + ylab('Standard Deviation')
g = g + theme(axis.text=element_text(size=16))
g = g + theme(axis.title=element_text(size=16,face='bold'))
g

ggsave('plots/sdsvsmean.png')

# plot with all_filter_rep is saved as sdsvsmean_unnormalized.png

# three cases 
a4b7_ss_ave = rowMeans(all_fn_rep[,all_filter.types %in% '1 a4b7+'])
a4b7_b_ave = rowMeans(all_fn_rep[,all_filter.types %in% '10 a4b7+'])
a4b7_b_sd = apply(all_fn_rep[,all_filter.types %in% '10 a4b7+'],1,sd)

ilcp_ss_ave = rowMeans(all_fn_rep[,all_filter.types %in% '1 ILCP'])
ilcp_b_ave = rowMeans(all_fn_rep[,all_filter.types %in% '10 ILCP'])
ilcp_b_sd = apply(all_fn_rep[,all_filter.types %in% '10 ILCP'],1,sd)

ltip_ss_ave = rowMeans(all_fn_rep[,all_filter.types %in% '1 LTiP'])
ltip_b_ave = rowMeans(all_fn_rep[,all_filter.types %in% '10 LTiP'])
ltip_b_sd = apply(all_fn_rep[,all_filter.types %in% '10 LTiP'],1,sd)

a4b7.df = data.frame(ss_ave=a4b7_ss_ave,b_ave=a4b7_b_ave,b_sd=a4b7_b_sd,
           symbol=rownames(all_fn_rep),type='a4b7+')

ilcp.df = data.frame(ss_ave=ilcp_ss_ave,b_ave=ilcp_b_ave,b_sd=ilcp_b_sd,
           symbol=rownames(all_fn_rep),type='ILCP')

ltip.df = data.frame(ss_ave=ltip_ss_ave,b_ave=ltip_b_ave,b_sd=ltip_b_sd,
           symbol=rownames(all_fn_rep),type='LTiP')

ssb.df = rbind(a4b7.df,ilcp.df,ltip.df)

# plot recapitulation
g = ggplot(ssb.df,aes(x=ss_ave,y=b_ave,color=type))
g = g + geom_point() + geom_errorbar(aes(ymin=b_ave-b_sd,ymax=b_ave+b_sd),width=1)
g = g + xlab('Single Cell Mean') + ylab('10 Cell Mean')
g = g + theme(axis.text=element_text(size=16))
g = g + theme(axis.title=element_text(size=16,face='bold'))
g

ggsave('plots/single_cell_vs_bulk.png')

# plots
# some type helpers
agg_types = all_filter.types
agg_types[all_filter.types %in% c("1 a4b7+","10 a4b7+")] = 'a4b7+'
agg_types[all_filter.types %in% c("1 ILCP","10 ILCP")] = 'ILCP'
agg_types[all_filter.types %in% c("1 LTiP","10 LTiP")] = 'LTiP'

bulkss_types = all_filter.types
bulkss_types[all_filter.types %in% c("1 a4b7+","1 ILCP","1 LTiP")] = 'Single Cell'
bulkss_types[all_filter.types %in% c("10 a4b7+","10 ILCP","10 LTiP")] = '10 Cell'

ss_pca = prcomp(t(all_fn_rep[gene_sds > 0,]),scale=T)

ss_pca.df = data.frame(ss_pca$x[,1:10])
ss_pca.df$NCells = bulkss_types
ss_pca.df$Type = agg_types

var_explained = (ss_pca$sdev)^2 / sum(ss_pca$sdev^2)

g = ggplot(ss_pca.df, aes(x=PC2,y=PC3,color=Type,shape=NCells)) + geom_point(size=4)
g = g + xlab('PC2 (6.89%)') + ylab('PC3 (4.02%)')
g = g + theme(axis.text=element_text(size=16))
g = g + theme(axis.title=element_text(size=16,face='bold'))
g

ggsave('plots/ss_pca.png')

## PCA just among single cells
ss_gene_sds = apply(all_fn_rep[,bulkss_types=='Single Cell'],1,sd)
all_fn_rep_var = all_fn_rep[ss_gene_sds > 0,]
ss_single_pca = prcomp(t(all_fn_rep_var[,bulkss_types=='Single Cell']),scale=T)

ss_single_pca.df = data.frame(ss_single_pca$x[,1:10])
ss_single_pca.df$Type = agg_types[bulkss_types=='Single Cell']

var_explained = (ss_single_pca$sdev)^2 / sum(ss_single_pca$sdev^2)

g = ggplot(ss_single_pca.df, aes(x=PC1,y=PC3,color=Type)) + geom_point(size=4)
g = g + xlab('PC1 (12.19%)') + ylab('PC3 (7.21%)')
g = g + theme(axis.text=element_text(size=16))
g = g + theme(axis.title=element_text(size=16,face='bold'))
g

ggsave('plots/ss_single_pca.png')

g = ggplot(ss_single_pca.df, aes(x=PC1,y=PC2,color=Type)) + geom_point(size=4)
g = g + xlab('PC1 (12.19%)') + ylab('PC2 (7.64%)')
g = g + theme(axis.text=element_text(size=16))
g = g + theme(axis.title=element_text(size=16,face='bold'))
g

ggsave('plots/ss_single_pca_PC12.png')

# sort(abs(ss_single_pca$rotation[,1]))

# heatmap
require(gplots)

png('plots/ss_heatmap.png',width=1500,height=1500)
heatmap.2(as.matrix(all_fn_rep),trace='none',density.info='none',scale='none',
          col=redgreen(75),srtCol=90,cexRow=1.5,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

png('plots/ss_heatmap_scale.png',width=1500,height=1500)
heatmap.2(as.matrix(all_fn_rep),trace='none',density.info='none',scale='row',symbreaks=F,
          col=redgreen(100),breaks=seq(-3,3,length.out=101),
          srtCol=90,cexRow=1.5,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# plotting correlation between cells
require(RColorBrewer)
bluescale = colorRampPalette(brewer.pal(9,"Blues"))(100)
png('plots/ss_cor_heatmap.png',width=1000,height=1000)
heatmap.2(cor(as.matrix(all_fn_rep)),Colv='Rowv',trace='none',density.info='none',col=bluescale,
          srtCol=90,cexRow=1.5,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# plot correlation between transcription factors
zero_var = apply(all_fn_rep,1,var) == 0
all_fn_rep_var = all_fn_rep[!zero_var,]
png('plots/gg_cor_heatmap.png',width=1200,height=1200)
heatmap.2(cor(as.matrix(t(all_fn_rep_var))),Colv='Rowv',trace='none',density.info='none',col=bluescale,
          srtCol=90,cexRow=1.5,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# sorting
# tf_cor = cor(as.matrix(t(all_fn_rep_var)))
# sort(tf_cor[rownames(tf_cor)=='Zbtb16',])
