# -------------------------------
# -------------------------------
# Code to run analyses in the corticostriatal paper
# 
# Gene expression links functionally coupled aspects of cortex and striatum
# Anderson, K.M., Krienen, F.M., Choi, E.Y.3, Reinen, J.M., Yeo, B.T., Holmes. A.J.
#
#
# Written by: Kevin M. Anderson
# Contact:    kevin.anderson@yale.edu 
# -------------------------------
# -------------------------------

library(data.table) 
library(WGCNA) 
library(psych) 
library(doBy) 
library(limma) 
library(plyr) 
library(gplots) 
library(ggplot2) 
library(biomartr) 
library(pbapply) 
library(biomaRt) 
library(edgeR) 
library(readxl) 
library(GEOquery) 


# Modify these filepaths for your local directory structure
# ----------------
afni.dir <- '/Users/kevinanderson/abin/'  # enter the abin directory for local install of AFNI
base.dir <- '/Users/kevinanderson/PHD/PROJECTS/2017_CORTICOSTRIATAL_NATCOMM/'  # enter the directory containing the script repository

# SET BASE DIRECTORY; source function library
# ----------------
function.lib <- paste(base.dir, 'scripts/function_library.R', sep = '')
source(function.lib)
# ----------------

# READ AHBA DATA
# ----------------
# (1) Microarray expression, 
# (2) Probes/sample information. 
# ----------------
load(paste(base.dir, 'data/AHBA/all_data.Rdata', sep = ''))
# ----------------

# SUBJECT LIST
# ----------------
donor.nums <- c('9861', '10021', '12876', '14380', '15496', '15697')
# ----------------

# Get info/sample data about all cortical and striatal samples, 
# store them in separate frames in the 'all_data' data structure
# ----------------
# CORTEX
cort_in     <- read.csv(paste(base.dir, 'reference_files/cort_regions.csv', sep = ''), header = FALSE) # Ontology IDs corresponding to cortical samples
cortex      <- as.numeric(cort_in$V1)
all_data    <- get_region_info(all_data=all_data, filenames=donor.nums, reg_IDs=cortex, name='cort')

# STRIATUM
striat_in   <- read.csv(paste(base.dir, 'reference_files/striat_regions.csv', sep = ''), header = FALSE) # Ontology  IDs corresponding to striatal samples
striatum    <- as.numeric(striat_in$V1)
all_data    <- get_region_info(all_data=all_data, filenames=donor.nums, reg_IDs=striatum, name='striat')
# ----------------

# pearson correlations between each 7-network striatal data to each 17-network cortical parcel
# ----------------
striat2cort.fcmri <- read.csv(paste0(base.dir, 'reference_files/choi7_yeo17_func.csv'), header=T)
# ----------------
# Read in the 114x114 cortical rs-fcMRI correlation matrix
# ----------------
cort2cort.fcmri <- read.csv(paste0(base.dir, 'reference_files/func_cor_mat.csv'), row.names = 1, header = FALSE)
colnames(cort2cort.fcmri) <- rownames(cort2cort.fcmri)
# ----------------

# Read previously defined info about the overlap of each sample to the cortical/striatal atlases
# ----------------
all_data <- readAtlasOverlap(all_data, filenames=donor.nums, atlas.dir=paste0(base.dir, 'atlas_overlap'))
# ----------------

# ----------------
# Read info about split label names for cortex, IDs, and color labels
# this splits the Yeo atlas into 57/114 spatially contiguous regions, depending on whether it's the 7- or 17-network atlas
# ----------------
atlas.key.17 <- read.table(paste0(base.dir, 'data/Yeo_JNeurophysiol11_SplitLabels/MNI152/17Networks_ColorLUT_freeview.txt'), col.names=c('ID','Name','R','G','B','A'))
atlas.key.7  <- read.table(paste0(base.dir, 'data/Yeo_JNeurophysiol11_SplitLabels/MNI152/7Networks_ColorLUT_freeview.txt'), col.names=c('ID','Name','R','G','B','A'))
# ----------------

# ----------------
# Naming information about the Choi striatal regions (e.g. Default Mode=7)
# ----------------
choi.names <- NULL
choi.names[['sev']]     <- as.character(read.csv(paste0(base.dir, 'reference_files/choi7_names.csv'), header=F)$V1)
choi.names[['sevteen']] <- as.character(read.csv(paste0(base.dir, 'reference_files/choi17_names.csv'), header=F)$V1)
# ----------------

# Find cortical parcels that contain a sample in the at least 2 of 6 subjects
# ----------------
# Across both hemispheres, all subjects, 7 networks.
use.regions.7 <- getCortRegions(all_data, donor.nums, atlas.num='7', thresh=2)
# ----------------

# Seperately for each donor, average expression of samples in the same cortical/striatal functional parcel
# ----------------
all_data <- avg_parcel_expression(all_data, 
                                  striatal.atlas='ChoiMNI152_7', # We opt for the 7 network parcellation in striatum to minimize functional regions with sparse sampling. 
                                  striatal.num='sev',
                                  cortical.atlas='splitLabel_7',
                                  cortical.num=51,
                                  cort.use.regions=use.regions.7)


# Array to match each of the 51 spatially contiguous parcels to a Yeo cortical network name (e.g. Parcel 23=Default)
# ------------------------------------
reg.names.27 <- as.character(atlas.key.7$Name[2:27]) # Left hemispheres
reg.names.51 <- as.character(atlas.key.7$Name[2:52]) # both hemispheres
reg2yeo.27   <- matrix(NA,26)
reg2yeo.51   <- matrix(NA,51)
for(choi in choi.names[['sev']]){
  cur.idxs     <- grep(choi, reg.names.27)
  cur.idxs.51 <- grep(choi, reg.names.51)
  # Treat TempPar as Default, like in original Yeo paper
  if (choi == 'Default'){
    cur.idxs     <- c(cur.idxs, grep('TempPar', reg.names.27))
    cur.idxs.51 <- c(cur.idxs.51, grep('TempPar', reg.names.51))
  }
  reg2yeo.27[cur.idxs] <- choi
  reg2yeo.51[cur.idxs.51] <- choi
}


# genes/regions in cortical analysis
# ----------------
gene.list    <- rownames(all_data[[1]]$all_cort_micros)
region.names <- c('Default', 'Cont', 'Limbic', 'VentAttn', 'DorsAttn', 'SomMot', 'Vis')
# ----------------


# Calculate differential expression of genes in 4 lh-hemisphere subjects 
# ------------------------------------------------------
lh.use.regions.7 <- getCortRegions(all_data, donor.nums[3:6], atlas.num='7', thresh=2)
diff.expr.list   <- ahba_diff_expression(all_data, 
                                         use.donors=donor.nums[3:6], 
                                         use.regions=lh.use.regions.7,
                                         reg2yeo=reg2yeo.27,
                                         dat2avg='cort_expr_nonorm',
                                         rest.networks=region.names)
cort.foldchange.n4 <- diff.expr.list[[1]]
n4.sig.genes       <- diff.expr.list[[2]]
# Save output to ./data directory
# ------------------------
n4.sig.genes <- unique(n4.sig.genes)
name.out     <- paste(base.dir, 'output_files/SUPPDATA_sheet2_n4_sig_cortical_genes.csv', sep = '')
write.table(x=n4.sig.genes, file=name.out, quote=TRUE, col.names=FALSE, row.names=FALSE)


# SUPPLEMENTAL DATA TABLE
# format the differential expression information into one big table
# --------------------------
out.matrix  <- NULL
col.headers <- NULL
for ( reg in region.names ) {
  out.dat     <- cbind(cort.foldchange.n4$stats[[reg]]$logFC, cort.foldchange.n4$stats[[reg]]$P.Value, cort.foldchange.n4$stats[[reg]]$adj.P.Val)
  out.matrix  <- cbind(out.matrix, out.dat)
  col.headers <- c(col.headers, c(paste(reg, 'logFC', sep='_'), paste(reg, 'pval', sep='_'), paste(reg, 'adj.pval', sep='_')))
}
colnames(out.matrix) <- col.headers
rownames(out.matrix) <- rownames(cort.foldchange.n4$stats[[reg]])
write.csv(x=out.matrix, file=paste0(base.dir, 'output_files/SUPPDATA_sheet1_n4_ahba_net7_cortical_fold_changes.csv'))

# Figure 2C
# Order of cortical parcels for plotting in a correlation grid. 
# 17-network parcels are used To increase resolution of cortico-cortical relationships. 
# --------------------------
twohemi.net17.regions <- getCortRegions(all_data, donor.nums[1:2], '17', 2)
plot.order    <- read.csv(paste0(base.dir, 'reference_files/17_network_plot_order_v2.txt'), header = FALSE)
plot.order$V1 <- rev(plot.order$V1)
all.idxs      <- NULL
reg.names     <- as.character(atlas.key.17$Name[2:dim(atlas.key.17)[1]])
for (reg in twohemi.net17.regions){
  cur.name <- reg.names[reg]
  cur.idx  <- which(plot.order$V1 %in% cur.name)
  all.idxs <- c(all.idxs, cur.idx)
}
ordered.use.regions <- twohemi.net17.regions[order(all.idxs)]
# --------------------------


# Average expression within each parcel that contains data from each bi-hemispheric subject. 
# ----------------
lh.use.regions.17 <- getCortRegions(all_data, donor.nums[1:2], atlas.num='17', thresh=2)
for ( donor in donor.nums[1:2]){
  print(paste('Averaging cort and striat data for: ', donor, sep = ''))
  
  # Set the modifiable atlas parameters, for instance you can run striatal 17/cortical 17
  # ----------------------
  cortical.atlas <- 'splitLabel_17'
  cortical.num   <- 114
  cort.use.regions <- lh.use.regions.17
  # Mean Normalized cortical expression values
  all_data[[donor]]$cort_expr_17         <- averageCortExpr(cort.use.regions, all_data[[donor]], cortical.num, 
                                                            all_data[[donor]][[cortical.atlas]], 'cort_meanNorm')
  # non-Mean Normalized cortical expression values
  all_data[[donor]]$cort_expr_nonorm_17  <- averageCortExpr(cort.use.regions, all_data[[donor]], cortical.num, 
                                                            all_data[[donor]][[cortical.atlas]], 'all_cort_micros')
}
# Network-associated genes used for the cortico-cortical correlations
# --------------------------
n4.sig.genes <- unique(n4.sig.genes)


# Make the mRNA cortico-cortical correlation matrix
# --------------------------
gene.idxs <- which(rownames(all_data[['10021']]$cort_expr) %in% unique(n4.sig.genes))
out.mats  <- NULL
# Calculate each gene co-experssion matrix seperately for each subject
for ( donor in donor.nums[1:2] ){
  cur.mat           <- all_data[[donor]]$cort_expr_17[gene.idxs, ordered.use.regions]
  untransformed     <- cor(cur.mat, method = 'spearman')
  #untransformed     <- cor(cur.mat, method = 'pearson')
  
  scaled  <- untransformed # z-tranfsorm spearman correlations
  
  z.corrs <- fisherz(untransformed[upper.tri(untransformed)]) # Upper portion of corr grid
  scaled[upper.tri(untransformed)] <- z.corrs
  
  z.corrs <- fisherz(untransformed[lower.tri(untransformed)])  # Lower portion of corr grid
  scaled[lower.tri(untransformed)] <- z.corrs
  
  # Label the rows/columns
  out.mats[[donor]] <- scaled
  rownames(out.mats[[donor]]) <- reg.names[ordered.use.regions]
  colnames(out.mats[[donor]]) <- reg.names[ordered.use.regions]
}
# Average the correlation matrices
avg.mat   <- (out.mats[[1]] + out.mats[[2]])/2 # average the z-transformed correlations of each donors
reg.names <- as.character(atlas.key.17$Name[2:length(atlas.key.17$Name)]) # add region names
rownames(avg.mat) <- reg.names[ordered.use.regions]
colnames(avg.mat) <- reg.names[ordered.use.regions]
# --------------------------


# Read  blue-black-red colormap
# --------------------------
mycol       <- rgb(read.table(paste0(base.dir, 'reference_files/rb_colormap.txt'), header = FALSE))
mrna.colors <- c(seq(-.6,.6,length=length(mycol)+1))
# --------------------------
# Plot the cortico-cortical mrna coexpression matrix
# --------------------------
diag(avg.mat) <- 0 # set diagonal to 0, for aesthetics
heatmap.2(as.matrix(avg.mat), trace ='none', col = mycol, Colv=FALSE, Rowv = FALSE, dendrogram = 'none',
          breaks=mrna.colors, symm=F,symkey=F,symbreaks=T, scale="none")
diag(avg.mat) <- 1
write.csv(file = paste(base.dir, 'output_files/corticocortical_mrna_correlations.csv', sep=''), x = avg.mat)
# --------------------------


# Figure 2B
# Plot the cortico-cortical rs-fcMRI matrix
# --------------------------
# reverse the rs-fcmri matrix, and flip horizontaly
mtx.tmp.v <- apply(cort2cort.fcmri, 2, rev)
cort2cort.fcmri <- apply(mtx.tmp.v, 1, rev)
# Select the rows that correspond to the mRNA data
# --------------------------
use_idxs        <- colnames(cort2cort.fcmri) %in% rownames(avg.mat)
plot.func       <- cort2cort.fcmri[use_idxs, use_idxs]
diag(cort2cort.fcmri) <- 0 
fc.colors <- c(seq(-.6,.6,length=length(mycol)+1))
heatmap.2(as.matrix(plot.func), trace ='none',col = mycol, Colv=FALSE, Rowv = FALSE, dendrogram = 'none',
          breaks=fc.colors, symm=F, symkey=F, symbreaks=T, scale="none")
diag(plot.func) <- 1
write.csv(file = paste(base.dir, 'output_files/corticocortical_rsfMRI_correlations.csv', sep=''), x = plot.func)
# --------------------------

# Produce a .csv file that will be converted to a paint file for plotting in Caret. 
# Figure 2A
# ------------------------
twohemi.use.regions <- getCortRegions(all_data, donor.nums[1:2], atlas.num = '17', thresh = 2)
bihemi.regions      <- matrix(-1,114,1)
bihemi.regions[twohemi.use.regions] <- 1 # usable regions are 1, everything else is -1
rownames(bihemi.regions) <- atlas.key.17$Name[2:115]
write.csv(file = paste(base.dir, '/output_files/bihemi_regions_59parcels.csv', sep = ''), x=bihemi.regions)
# ------------------------


# Global correspondence of mrna/fcmri
# ----------------
func.1d.arr <- plot.func[upper.tri(plot.func)]
mrna.1d.arr <- avg.mat[upper.tri(avg.mat)]
cor.test(func.1d.arr, mrna.1d.arr)
# ----------------


# Figure 2C
# ---------------------
# Determine if wihtin-network genetic correlations are higher than between network...
tmp.reg.names <- c('Default|TempPar', "Cont", "Limbic", "VentAttn", "DorsAttn", "SomMot", "Vis")
cor.table.out <- NULL
for ( donor in donor.nums[1:2]) {
  out.mat <- out.mats[[donor]]
  
  for (reg in tmp.reg.names) {
    # grab the within-network correlations
    # -------------------
    cur.idxs   <- grep(reg, rownames(out.mat))
    cur.within <- out.mat[cur.idxs, cur.idxs][upper.tri(out.mat[cur.idxs, cur.idxs])]
    avg.within <- mean(cur.within)
    names(avg.within) <- reg
    # get the average corr each of the non-corresponding networks
    o.regs    <- tmp.reg.names[which(!tmp.reg.names %in% reg)]
    out.corrs <- NULL
    for (oreg in o.regs){
      cur.o.idxs <- grep(oreg, rownames(out.mat))
      cur.o.avg  <- mean(colMeans(out.mat[cur.idxs, cur.o.idxs]))
      out.corrs  <- c(out.corrs, cur.o.avg)
    }
    names(out.corrs) <- o.regs
    cors.out <- c(avg.within, out.corrs)
    
    # Organize the output
    # ---------------
    cur.reg   <- rep(reg, 1, 7)
    cur.donor <- rep(donor, 1, 7)
    groups    <- c('within', rep('between', 1,6))
    cur.out   <- cbind(cur.reg, cur.donor, cors.out, groups)
    cor.table.out <- rbind(cor.table.out, cur.out)
  }
}
rownames(cor.table.out) <- NULL
cor.table.out           <- as.data.frame(cor.table.out)
colnames(cor.table.out) <- c('region', 'donor', 'corr', 'group')
cor.table.out$corr      <- as.numeric(as.character(cor.table.out$corr))
cor.table.out$group  <- factor(as.character(cor.table.out$group), levels=c('within', 'between'))

# Get data ready for plotting
# --------------------
dat.df <- NULL
for ( reg in tmp.reg.names ){
  within.idxs <- intersect(which(cor.table.out$region == reg), which(cor.table.out$group == 'within'))
  within.mean <- mean(cor.table.out[within.idxs, ]$corr)
  within.se   <- sd(cor.table.out[within.idxs, ]$corr)/sqrt(2)
  
  # donor 1
  between.idxs.d1 <- intersect(intersect(which(cor.table.out$region == reg), which(cor.table.out$group == 'between')), which(cor.table.out$donor == '9861'))
  between.mean.d1 <- mean(cor.table.out[between.idxs.d1, ]$corr)
  between.sd.d1   <- sd(cor.table.out[between.idxs.d1, ]$corr)
  # donor 2
  between.idxs.d2 <- intersect(intersect(which(cor.table.out$region == reg), which(cor.table.out$group == 'between')), which(cor.table.out$donor == '10021'))
  between.mean.d2 <- mean(cor.table.out[between.idxs.d2, ]$corr)
  between.sd.d2   <- sd(cor.table.out[between.idxs.d2, ]$corr)  
  
  between.mean <- mean(c(between.mean.d1, between.mean.d2))
  between.se   <- sd(c(between.mean.d1, between.mean.d2))/sqrt(2)
  
  cur.row <- c(within.mean, within.se, between.mean, between.se)
  dat.df <- rbind(dat.df, cur.row)
}
colnames(dat.df) <- c('wmean', 'wse', 'bmean', 'bse')
rownames(dat.df) <- tmp.reg.names
dat.df <- as.data.frame(dat.df)


# Plot network-specific cortical expression in the 2 bi-hemispheric subjects
# --------------------------
dat.df  <- dat.df[c('Default|TempPar','Cont','Limbic','VentAttn','DorsAttn','SomMot','Vis'),]
plot.df <- as.data.frame(cbind(c(dat.df$wmean, dat.df$bmean), c(dat.df$wse, dat.df$bse), 
                               c(rownames(dat.df), rownames(dat.df)), c(rep('awithin',7), rep('zbetween',7))))
colnames(plot.df) <- c('expr','se','net', 'group')
plot.df$expr      <- as.numeric(as.character(plot.df$expr))
plot.df$se        <- as.numeric(as.character(plot.df$se))
plot.df$net       <- factor(plot.df$net, levels=c('Default|TempPar','Cont','Limbic','VentAttn','DorsAttn','SomMot','Vis'))
#plot_df$expr[c(1,8)] <- 0
#plot_df$se[c(1,8)]  <- 0
p <- ggplot(plot.df, aes(fill=group, y=expr, x=net))
limits <- aes(ymax = expr + se, ymin=expr - se)
dodge  <- position_dodge(width=0.9)
dot.df <- summaryBy(corr ~ region + group + donor, cor.table.out)
ggplot(data=plot.df, aes(fill=group, y=expr, x=net)) + geom_bar(position=dodge, stat="identity") + 
  geom_errorbar(limits, width=0.25) + ylim(c(-.5, 1)) +
  geom_dotplot(data=dot.df, aes(x=region, y=corr.mean, fill=group), binaxis='y', dotsize = .4, stackdir='center', position=position_dodge(0.8))
# --------------------------


# Get Average cortical expression for each subject, in 6 subjects. 
# --------------------------
# 'cort_expr' is used for the plotting
region.names         <- c('Vis','VentAttn','DorsAttn','SomMot','Limbic','Cont','Default')
use.n6.net7.regions  <- getCortRegions(all_data, donor.nums, '7', 2)
plot.out             <- averageWithinCortNetworks(all_data, use.n6.net7.regions, reg2yeo.51, donor.nums, 'cort_expr', region.names)
plot.n6.cort.expr    <- plot.out[[1]]
plot.n6.cort.regions <- plot.out[[2]]
# Plot the expression of exemplar genes across cortical regions -- used for supplementary figures
# --------------------------
makeGenePlot(all_data, 'SSTR1', gene.list, plot.n6.cort.expr, plot.n6.cort.regions, donor.nums, ylimit=2)
makeGenePlot(all_data, 'PVALB', gene.list, plot.n6.cort.expr, plot.n6.cort.regions, donor.nums, ylimit=2)


# --------------------------
# Calculate cortical differential expression across all 6 subjects. 
# --------------------------
region.names         <- c('Vis','VentAttn','DorsAttn','SomMot','Limbic','Cont','Default')
use.n6.net7.regions  <- getCortRegions(all_data, donor.nums, '7', 2)
diff.expr.list <- ahba_diff_expression(all_data, 
                                       use.donors=donor.nums, 
                                       use.regions=use.n6.net7.regions,
                                       reg2yeo=reg2yeo.51,
                                       dat2avg='cort_expr_nonorm',
                                       rest.networks=region.names)
cort.foldchanges.n6 <- diff.expr.list[[1]]
sig.n6.genes        <- diff.expr.list[[2]]

# SUPPLEMENTAL DATA TABLE 
# --------------------------
out.matrix  <- NULL
col.headers <- NULL
for ( reg in region.names ) {
  out.dat     <- cbind(cort.foldchanges.n6$stats[[reg]]$logFC, cort.foldchanges.n6$stats[[reg]]$P.Value, cort.foldchanges.n6$stats[[reg]]$adj.P.Val)
  out.matrix  <- cbind(out.matrix, out.dat)
  col.headers <- c(col.headers, c(paste(reg, 'logFC', sep='_'), paste(reg, 'pval', sep='_'), paste(reg, 'adj.pval', sep='_')))
}
colnames(out.matrix) <- col.headers
rownames(out.matrix) <- gene.list
name.out <- paste(base.dir, 'output_files/SUPPDATA_sheet3_n6_ahba_net17_cortical_fold_changes.csv', sep = '')
write.csv(x = out.matrix, file = name.out)
# -------------------------------


# Genes that are significantly differentially expressed in one of the cortical networks, n=6
# ------------------------------------
write.me    <- NULL
write.genes <- vector()
write.nums  <- NULL
for (choi in region.names){
  print(choi)
  genes.tmp   <- cort.foldchanges.n6[['q01_genes']][[choi]]
  genes.in    <- paste(genes.tmp, collapse = ',')
  write.genes <- rbind.fill.matrix(write.genes, t(genes.tmp))
  write.nums  <- c(write.nums, length(genes.tmp))
  write.me    <- rbind(write.me, c(choi, length(genes.tmp), paste(genes.tmp, collapse = ',')))
}
# Write the significant genes in readable format
write.genes[is.na(write.genes)] <- ''
t.write.genes <- t(write.genes)
write.me      <- rbind(region.names, rbind(write.nums, t.write.genes))
write.csv(x=write.me, file=paste0(base.dir, 'output_files/SUPPDATA_sheet5_n6_cortical_fdr_sig_genes.csv'), row.names=FALSE)

# Write the genes that passed the sig threshold
sig.n6.genes <- unique(sig.n6.genes)
write.table(x=sig.n6.genes, file=paste0(base.dir, 'output_files/SUPPDATA_sheet4_n6_ahba_pos_fdr01_genes.csv'), row.names=FALSE, col.names=FALSE, quote=TRUE)
# ---------------------


# Plot striatal expression
# ---------------------
# Valid striatal networks (i.e. data exists)
use.striat.networks <- c('Default','Cont','Limbic','VentAttn','SomMot')
plot.striat.out     <- averageWithinStriatNetworks(all_data, 
                                                   donor.nums, 
                                                   use.striat.networks, 
                                                   type='striat_meanNorm', 
                                                   net_names=choi.names$sev, 
                                                   atlas_field='ChoiMNI152_7')
plot.striat.expr    <- plot.striat.out[[1]]
plot.striat.regions <- plot.striat.out[[2]]
plot.donor.array    <- plot.striat.out[[3]]
# Plot the expression of exemplar genes across cortical regions -- used for supplementary figures
# --------------------------
use.striat.networks <- c('Limbic','Default','Cont','VentAttn','SomMot')
ahbaPlotStriatalExpr(plot.striat.expr, plot.striat.regions, plot.donor.array, use.striat.networks, 'SSTR1', ylimit = 3)
ahbaPlotStriatalExpr(plot.striat.expr, plot.striat.regions, plot.donor.array, use.striat.networks, 'PVALB', ylimit = 4)


# Identify genes that are differentially expressed within each striatal network
# ---------------------
# Get Average striatal expression for each subject, in 6 subjects. 
# --------------
use.striat.networks <- c('Default','Cont','Limbic','VentAttn','SomMot')
striat.out     <- averageWithinStriatNetworks(all_data, donor.nums, use.striat.networks, 'all_striat_micros', choi.names$sev, 'ChoiMNI152_7')
striat.expr    <- striat.out[[1]]
striat.regions <- striat.out[[2]]
donor.array    <- striat.out[[3]]
# Construct design matrix
# ------------------------------
striat.design           <- model.matrix(~0+as.factor(striat.regions))
colnames(striat.design) <- gsub('as.factor\\(striat.regions\\)', '', colnames(striat.design))
corfit <- duplicateCorrelation(striat.expr, striat.design, block=donor.array)
# Calculate differential expression for each striatal subregion
# ------------------------------
striat.foldchanges.n6 <- NULL
write.nums  <- NULL
write.genes <- vector()
for ( choi in use.striat.networks ) {
  print(choi)
  mult.term    <- round(1/(length(use.striat.networks)-1),6)
  o.chois      <- use.striat.networks[use.striat.networks != choi]
  cur.contrast <- paste('1*', choi, '-', mult.term, '*', paste(o.chois, collapse = paste('-', mult.term, '*', sep = '')), sep = '')
  cmtx         <- makeContrasts(cur.contrast, levels=striat.design ) # Make a contrast matrix
  tmplm        <- lmFit( striat.expr, striat.design, block=donor.array, correlation=corfit$consensus.correlation) # Fit the model to the data
  fit          <- eBayes( contrasts.fit( tmplm, cmtx ) )
  striat.foldchanges.n6[['fit_df']][[choi]] <- fit
  tmp          <- topTable(fit, number=Inf)
  striat.foldchanges.n6[['stats']][[choi]] <- tmp[order(rownames(tmp)),]
  pos_coefs    <- which(striat.foldchanges.n6[['stats']][[choi]]$logFC > 0)
  q.values     <- which(striat.foldchanges.n6[['stats']][[choi]]$adj.P.Val <= .01)
  genes_tmp    <- rownames(striat.foldchanges.n6[['stats']][[choi]])[intersect(q.values, pos_coefs)]
  print(length(genes_tmp))
  striat.foldchanges.n6[['q01_genes']][[choi]] <- genes_tmp
  
  write.genes  <- rbind.fill.matrix(write.genes, t(genes_tmp))
  write.nums   <- c(write.nums, length(genes_tmp))
}
write.csv(x=striat.foldchanges.n6$fold_changes, file=paste0(base.dir, 'output_files/n6_ahba_striatal_fold_changes.csv'))
# ------------------------------------------------------


# Identify genes enriched in both the striatal and cortical aspects of the limbic loop
# ------------------------------------
limbic.genes  <- intersect(striat.foldchanges.n6$q01_genes$Limbic, cort.foldchanges.n6$q01_genes$Limbic)
limbic.idxs   <- which(all_data$`9861`$probes_collapse$gene_symbol %in% limbic.genes)
limbic.entrez <- all_data$`9861`$probes_collapse$entrez_id[limbic.idxs]
# Somato/motor loop
sommot.genes  <- intersect(striat.foldchanges.n6$q01_genes$SomMot, cort.foldchanges.n6$q01_genes$SomMot)
sommot.idxs   <- which(all_data$`9861`$probes_collapse$gene_symbol %in% sommot.genes)
sommot.entrez <- all_data$`9861`$probes_collapse$entrez_id[sommot.idxs]


# Determine the number of enriched genes that do are not annotated at all
# ------------------------------------
# LIMBIC
GO_tbl <- getGO(organism = "Homo sapiens", 
                genes    = limbic.entrez,
                filters  = "entrezgene")
length(limbic.entrez) - length(unique(GO_tbl$entrezgene))
limbic.unannotated_genes <- limbic.entrez[!limbic.entrez %in% unique(GO_tbl$entrezgene)]
# SOMATO/MOTOR
GO_sommot <- getGO(organism = "Homo sapiens", 
                   genes    = sommot.entrez,
                   filters  = "entrezgene")
length(sommot.entrez) - length(unique(GO_sommot$entrezgene))
sommot.unannotated_genes <- sommot.entrez[!sommot.entrez %in% unique(GO_sommot$entrezgene)]
# ------------------------------------


# Run hypergeometric tests - significant values indicate that a greater than expected number of genes overlap between cortex/striatum for a network
# Variable names are meant to mimick the popular example of drawing black/white balls from an urn. 
# ------------------------------------
hyper.sig.array <- NULL
write.genes <- vector()
write.nums  <- NULL
for (choi in use.striat.networks){
  sig.striatal.genes <- striat.foldchanges.n6$q01_genes[[choi]]
  sig.cortical.genes <- cort.foldchanges.n6$q01_genes[[choi]]
  
  num.overlap        <- length(intersect(sig.striatal.genes, sig.cortical.genes))
  num.genes          <- length(cort.foldchanges.n6$fit_df$Limbic$coefficients)
  white.genes        <- length(sig.striatal.genes)
  black.genes        <- num.genes-length(sig.striatal.genes)
  genes.drawn        <- length(sig.cortical.genes)
  print(choi)
  print(num.overlap)
  
  p.val <- phyper(num.overlap-1, white.genes, black.genes, genes.drawn, lower.tail = F)
  hyper.sig.array <- c(hyper.sig.array, p.val)
  
  overlap.genes <- intersect(sig.striatal.genes, sig.cortical.genes)
  if ( is.null(overlap.genes) ){
    overlap.genes <- ''
  }
  write.genes   <- rbind.fill.matrix(write.genes, t(overlap.genes))
  write.nums    <- c(write.nums, length(overlap.genes))
}
write.genes[is.na(write.genes)] <- ''
t.write.genes <- t(write.genes)
write.me  <- rbind(use.striat.networks, rbind(write.nums, t.write.genes))
base      <- paste(base.dir, 'output_files/SUPPDATA_sheet8_sig_cortexStriat_ahba_n6.csv', sep = '')
write.csv(x = write.me, file = base, row.names=FALSE)
# ------------------------------------

# SUPPLEMENTAL DATA TABLE 
# full information on striatal fold change estimates 
# --------------------------
out.matrix  <- NULL
col.headers <- NULL
for ( reg in use.striat.networks ) {
  out.dat     <- cbind(striat.foldchanges.n6$stats[[reg]]$logFC, striat.foldchanges.n6$stats[[reg]]$P.Value, striat.foldchanges.n6$stats[[reg]]$adj.P.Val)
  out.matrix  <- cbind(out.matrix, out.dat)
  col.headers <- c(col.headers, c(paste(reg, 'logFC', sep='_'), paste(reg, 'pval', sep='_'), paste(reg, 'adj.pval', sep='_')))
}
colnames(out.matrix) <- col.headers
rownames(out.matrix) <- gene.list
write.csv(x = out.matrix, file = paste0(base.dir, 'output_files/SUPPDATA_sheet6_n6_ahba_net7_striatal_fold_changes.csv'))

# Write the significant genes in a readable format
write.genes[is.na(write.genes)] <- ''
t.write.genes <- t(write.genes)
write.me <- rbind(use.striat.networks, rbind(write.nums, t.write.genes))
write.csv(x = write.me, file = paste0(base.dir, 'output_files/SUPPDATA_sheet7_n6_ahba_striatal_fdr01_genes.csv'), row.names=FALSE)
# ------------------------------------


# correlate gene-wise fold change from cortex/striatum
# -------------------
cor.test(cort.foldchanges.n6$fit_df$Limbic$coefficients, striat.foldchanges.n6$fit_df$Limbic$coefficients)
cor.test(cort.foldchanges.n6$fit_df$SomMot$coefficients, striat.foldchanges.n6$fit_df$SomMot$coefficients)
cor.test(cort.foldchanges.n6$fit_df$SomMot$coefficients, striat.foldchanges.n6$fit_df$Limbic$coefficients)
cor.test(cort.foldchanges.n6$fit_df$Limbic$coefficients, striat.foldchanges.n6$fit_df$SomMot$coefficients)

# Supplemental Figures -- genome-wide correlation of fold change scores
# -------------------
plotSuppFigure3(cort.foldchanges.n6, striat.foldchanges.n6)


# Get the within/between network expression of cortical network genes in striatum
# --------------------
striat.out.ind <- averageWithinStriatNetIndividual(all_data, 
                                                   donor.nums, 
                                                   use.striat.networks, 
                                                   type='striat_meanNorm', 
                                                   net_names=choi.names$sev, 
                                                   atlas_field='ChoiMNI152_7')
all.striat.out <- striatExprOfCortGeneSets(cort.foldchanges.n6, striat.out.ind)
aov.dat.striat <- all.striat.out[all.striat.out$network %in% c('Default','Cont','Limbic','VentAttn','SomMot'),]

aov.model  <- aov(expr ~ group + network + Error(donor/group), aov.dat.striat)
summary(aov.model)
summ.mat <- summaryBy(expr ~ group + network + Error(donor/group), aov.dat.striat)
within   <- NULL
between  <- NULL
for (don in levels(summ.mat$donor)){
  don.idxs <- which(summ.mat$donor == don)
  sub.data <- summ.mat[summ.mat$donor == don,]
  within   <- c(within, mean(sub.data$expr.mean[sub.data$group == 'within']))
  between  <- c(between, mean(sub.data$expr.mean[sub.data$group == 'between']))
}
mean(within)
sd(within)/sqrt(length(within))
mean(between)
sd(between)/sqrt(length(between))
# Post-hoc tests
# ------------------
for ( reg in use.striat.networks ){
  cat('*\n*\n*\n*\n*\n')
  print(reg)
  cur.dat     <- all.striat.out[all.striat.out$network == reg,]
  ttest.model <- aov(expr ~ group + Error(donor/group), cur.dat)
  print(summary(ttest.model))
}
# Test whether Som/Mot cortical genes are expressed more in VentAttn striatum
ventattn.aov.dat.striat <- all.striat.out[all.striat.out$network %in% c('VentAttn_SomMot'),]
ventattn.aov.model  <- aov(expr ~ group + Error(donor/group), ventattn.aov.dat.striat)
summary(ventattn.aov.model)


# FIGURE 3B - Plotting
# --------------------
plot.striat.out <- NULL
for (reg in c(use.striat.networks, 'VentAttn_SomMot')){
  cur.idxs <- intersect(which(all.striat.out$network == reg), which(all.striat.out$group == 'within'))
  cur.dat  <- all.striat.out[cur.idxs,]
  cur.within.mean <- mean(cur.dat$expr)
  cur.within.se   <- sd(cur.dat$expr)/sqrt(6)
  betweens <- NULL
  for (donor in donor.nums){
    reg.idxs   <- which(all.striat.out$network == reg)
    group.idxs <- which(all.striat.out$group == 'between')
    donor.idxs <- which(all.striat.out$donor == donor)
    betweens   <- c(betweens, mean(all.striat.out$expr[intersect(donor.idxs, intersect(reg.idxs, group.idxs))]))
  }
  cur.betw.mean   <- mean(betweens)
  cur.betw.se     <- sd(betweens)/sqrt(6)
  cur.row         <- c(cur.within.mean, cur.within.se, cur.betw.mean, cur.betw.se)
  plot.striat.out <- rbind(plot.striat.out, cur.row)
}
colnames(plot.striat.out) <- c('wmean','wse','bmean','bse')
rownames(plot.striat.out) <- c(use.striat.networks, 'VentAttn_SomMot')
dat.df <- plot.striat.out
# --------------------
# Plot striatal expression values for genes defined at cortex. 
dat.df <- as.data.frame(dat.df)
rownames(dat.df) <- c("Default","Cont","Limbic","VentAttn","SomMot",'VentAttn_SomMot')
plot_df <- as.data.frame(cbind(c(dat.df$wmean, dat.df$bmean), 
                               c(dat.df$wse, dat.df$bse), 
                               c(rownames(dat.df), rownames(dat.df)),
                               c(rep('within',6), rep('zbetween',6))))
colnames(plot_df) <- c('cor','se','net', 'group')
plot_df$cor       <- as.numeric(as.character(plot_df$cor))
plot_df$se        <- as.numeric(as.character(plot_df$se))
p <- ggplot(plot_df, aes(fill=group, y=cor, x=net))
plot_df$group <- as.factor(plot_df$group)
limits <- aes(ymax = cor + se, ymin=cor - se)
dodge  <- position_dodge(width=0.9)
plot_df$net <- factor(plot_df$net, levels = c("Default","Cont","Limbic","VentAttn","SomMot",'VentAttn_SomMot'))
dot.plot <- summaryBy(expr ~ + group + donor + network, all.striat.out)
dot.plot$group <- factor(dot.plot$group, levels=c('within','between'))
ggplot(data=plot_df, aes(fill=group, y=cor, x=net)) + 
  geom_bar(position=dodge, stat="identity") +
  scale_fill_manual(values=c('gray0','gray50', 'gray100')) +
  geom_errorbar(limits, width=0.25) + ylim(-.25, .5) +
  geom_dotplot(data=dot.plot, aes(x=network, y=expr.mean, fill=group), 
               binaxis='y', dotsize=.4, stackdir='center', position=position_dodge(0.8))
# --------------------


# Determine which cortical parcels are present across in at least 2 donor
# ---------------------------------------------------------
n17.use.regions  <- getCortRegions(all_data, donor.nums, '17', 2)
n17.use.regions  <- sort(n17.use.regions)
# -----------------
# Get the average expression of each 17-network parcel, across all six subjects. 
lh.use.regions.17 <- getCortRegions(all_data, donor.nums, atlas.num='17', thresh=2)
for ( donor in donor.nums){
  print(paste('Averaging cort and striat data for: ', donor, sep = ''))
  
  # Set the modifiable atlas parameters, for instance you can run striatal 17/cortical 17
  # ----------------------
  cortical.atlas <- 'splitLabel_17'
  cortical.num   <- 114
  cort.use.regions <- lh.use.regions.17
  # Mean Normalized cortical expression values
  all_data[[donor]]$cort_expr_17_73parcels <- averageCortExpr(cort.use.regions, all_data[[donor]], cortical.num, 
                                                              all_data[[donor]][[cortical.atlas]], 'cort_meanNorm')
}


# Genetic correlation between each cortical and striatal region, LH/RH, seperately for each donor
# -----------------
targ.idxs <- which(rownames(all_data[[donor]]$cort_expr) %in% sig.n6.genes)
for (donor in donor.nums){
  n.corts            <- 114 #
  n.striats          <- dim(all_data[[donor]]$striat_expr)[2] # this will vary by subj
  
  # use_regions refers to the cortical areas that have enough data from each subject
  corticostriat.cors           <- matrix(NA, n.corts, n.striats)
  colnames(corticostriat.cors) <- colnames(all_data[[donor]]$striat_expr)
  
  donor.regions <- getCortRegions(all_data, donor, atlas.num='17', thresh=1)
  cur.expr      <- averageCortExpr(cort.use.regions, all_data[[donor]], cortical.num, all_data[[donor]][[cortical.atlas]], 'cort_meanNorm')
  for ( cort in donor.regions ) {
    cur.cort <- cur.expr[targ.idxs, cort]
    if ( length(cur.cort[is.na(cur.cort)]) == 0 ) { # make sure the cortical area exists for this subject
      # Iterate over every striatal region
      for ( striat in 1:n.striats ){
        cur.striat                       <- all_data[[donor]]$striat_expr[targ.idxs, striat]
        corticostriat.cors[cort, striat] <- cor(cur.cort, cur.striat, method = 'spearman')
      }
      
    } else {
      for ( striat in 1:n.striats ){
        corticostriat.cors[cort, striat] <- NA
      }
    }
  }
  all_data[[donor]]$corticostriat.cors           <- as.data.frame(corticostriat.cors)
  rownames(all_data[[donor]]$corticostriat.cors) <- atlas.key.17$Name[2:115]
  rownames(corticostriat.cors) <- atlas.key.17$Name[2:115]
  
  # z-transform correlation values
  for (idx in 1:n.striats){
    corrs   <- corticostriat.cors[,idx]
    z.corrs <- fisherz(corrs[!is.na(corrs)])
    corrs[!is.na(corrs)] <- z.corrs
    corticostriat.cors[,idx] <- corrs
  }
  # output csv files with correlation values for plotting in caret
  # -----------------
  filename = paste(base.dir, 'caret_files/lhrh_', donor, '_n6_mrna_corrvals.csv', sep = '')
  write.csv(x=corticostriat.cors, file = filename)
}

# Average all of the cortico-striatal correlations across donors (after z-transform), for plotting
# and comparison with rs-fcMRI data
# -----------------------
avg.cort.striat.express <- NULL
min.striats <- c('Default','Cont','Limbic','VentAttn','SomMot')
for ( reg in min.striats ){
  expr.arr <- NULL
  cur.name <- reg
  for ( donor in donor.nums ){
    if (cur.name %in% colnames(all_data[[donor]]$corticostriat.cors) ){
      corrs   <- as.matrix(all_data[[donor]]$corticostriat.cors[[cur.name]])
      z.corrs <- fisherz(corrs[!is.na(corrs)])
      corrs[!is.na(corrs)] <- z.corrs
      expr.arr <- cbind(expr.arr, as.matrix(all_data[[donor]]$corticostriat.cors[[cur.name]]))
    }
  }
  avg.cort.striat.express <- cbind(avg.cort.striat.express, rowMeans(expr.arr, na.rm = TRUE))
}
# Add column names and convert to data.frame
colnames(avg.cort.striat.express) <- min.striats
avg.cort.striat.express.df        <- as.data.frame(avg.cort.striat.express)
rownames(avg.cort.striat.express.df) <- atlas.key.17$Name[2:115]
rownames(striat2cort.fcmri)    <- atlas.key.17$Name[2:115]
# Get rid of missing rows
avg.cort.striat.express.df.nonan  <- avg.cort.striat.express.df[which(!is.na(avg.cort.striat.express.df$Default)),]
striat2cort.fcmri.nonan        <- striat2cort.fcmri[which(!is.na(avg.cort.striat.express.df$Default)),]
rownames(striat2cort.fcmri.nonan) <- rownames(avg.cort.striat.express.df.nonan)
# Write output
# -----------------
filename = paste(base.dir, 'caret_files/lhrh_n6_mrna_corrvals.csv', sep = '')
write.csv(x=cbind(avg.cort.striat.express.df, 1:dim(avg.cort.striat.express.df)[1]), file = filename)
filename = paste(base.dir, 'caret_files/lhrh_n6_rsfc_corrvals.csv', sep = '')
write.csv(x = cbind(striat2cort.fcmri.nonan, 1:dim(striat2cort.fcmri.nonan)[1]), file = filename)
# -----------------


# Figure 2D/E CORTICO-CORTICAL Correlations
# Average cortical correlations for each bihemispheric cortical ROI - 'Limbic_OFC' and 'SomMotA'
# -----------------
# rs-fcMRI
# -----------------
export.order <- rownames(avg.cort.striat.express.df)
for (reg in c('Limbic_OFC', 'SomMotA')){
  export.cors  <- NULL
  region.cors  <- colMeans(cort2cort.fcmri[grep(reg, rownames(cort2cort.fcmri)),])
  for (exp in export.order){
    idx <- grep(paste('^', exp, '$', sep=''), names(region.cors))
    if ( length(idx) > 1 ){
      print(exp)
    }
    if (length(idx) == 0){
      export.cors <- c(export.cors, 'NA')
    } else {
      export.cors <- c(export.cors, as.character(region.cors[idx]))
    }
  }
  write.table(x=export.cors, file=paste0(base.dir, 'caret_files/rsfcmri_avg_', reg, '_cortical_corrs.csv'), sep=',', col.names = TRUE)
}
# -----------------
# Gene Expression
# -----------------
export.order   <- rownames(avg.cort.striat.express.df)
cort2cort.cors <- NULL
for (reg in c('Limbic_OFC', 'SomMotA')){
  region.cors  <- colMeans(avg.mat[grep(reg, rownames(avg.mat)),]) # avg.mat is the bi-hemi only cortico-cortical correlation matrix computed above
  export.cors  <- NULL
  for (exp in export.order){
    idx <- grep(paste('^', exp, '$', sep=''), names(region.cors))
    if ( length(idx) > 1 ){
      print(exp)
    }
    if (length(idx) == 0){
      export.cors <- c(export.cors, 'NA')
    } else {
      export.cors <- c(export.cors, as.character(region.cors[idx]))
    }
  }
  write.table(x=export.cors, file=paste0(base.dir, 'caret_files/mRNA_avg_', reg, '_cortical_corrs.csv'), sep=',', col.names = TRUE)
  cort2cort.cors[[reg]] <- export.cors
}
# -----------------


# Get the average correlation of each striatal network (column) to each cortical network (rows); 
# -----------------
cors.out <- matrix(NA, 7,5)
colnames(cors.out) <- c('Default','Cont','Limbic','VentAttn','SomMot')
rownames(cors.out) <- c('Default','Cont','Limbic','VentAttn','DorsAttn','SomMot','Vis')
colct <- 1
for (net in c('Default','Cont','Limbic','VentAttn','SomMot')){
  rowct <- 1
  cur.vals    <- avg.cort.striat.express.df.nonan[[net]]
  for (net2 in c('Default','Cont','Limbic','VentAttn','DorsAttn','SomMot','Vis')) {
    within.idxs <- grep(net2, rownames(avg.cort.striat.express.df.nonan))
    if (net2 == 'Default'){
      tmp <- grep('TempPar', rownames(avg.cort.striat.express.df.nonan))
      within.idxs <- c(within.idxs, tmp)
    }
    cur.mean.cor <- mean(cur.vals[within.idxs])
    cors.out[rowct, colct] <- cur.mean.cor
    rowct <- rowct + 1
  }
  colct <- colct + 1
}
cor.test(avg.cort.striat.express.df.nonan$Limbic, striat2cort.fcmri.nonan$Limbic)
cor.test(avg.cort.striat.express.df.nonan$Default, striat2cort.fcmri.nonan$Default)
cor.test(striat2cort.fcmri.nonan$Limbic, striat2cort.fcmri.nonan$Default)

cor.test(avg.cort.striat.express.df.nonan$SomMot, striat2cort.fcmri.nonan$SomMot)
cor.test(avg.cort.striat.express.df.nonan$Cont, striat2cort.fcmri.nonan$Cont)
cor.test(avg.cort.striat.express.df.nonan$VentAttn, striat2cort.fcmri.nonan$VentAttn)
cor.test(avg.cort.striat.express.df.nonan$SomMot, avg.cort.striat.express.df.nonan$VentAttn)
# -----------------
# Correlate the maps derived from limbic/sommot striatal seeds to the limbic/sommot cortical sseds
cort.limbicOFC  <- cort2cort.cors
use.idxs.cort   <- which(!cort2cort.cors$Limbic_OFC == 'NA')
use.idxs.striat <- which(!is.nan(avg.cort.striat.express.df$Limbic))
cor.test(avg.cort.striat.express.df$Limbic[intersect(use.idxs.cort, use.idxs.striat)], as.numeric(cort2cort.cors$Limbic_OFC[intersect(use.idxs.cort, use.idxs.striat)]))
cor.test(avg.cort.striat.express.df$SomMot[intersect(use.idxs.cort, use.idxs.striat)], as.numeric(cort2cort.cors$SomMotA[intersect(use.idxs.cort, use.idxs.striat)]))
# -----------------------
# Get average correlations of sommot cortex, from sommot striatal seed
mean(avg.cort.striat.express.df.nonan$SomMot[grep('SomMot', rownames(avg.cort.striat.express.df.nonan))])
sd(avg.cort.striat.express.df.nonan$SomMot[grep('SomMot', rownames(avg.cort.striat.express.df.nonan))]) / sqrt(length(avg.cort.striat.express.df.nonan$SomMot[grep('SomMot', rownames(avg.cort.striat.express.df.nonan))]))
# -----------------------


# Write files to be converted into caret surface images
# -----------------------
caret.dir <- paste(base.dir, 'caret_files', sep = '')
write.csv(file=paste0(caret.dir, '/SPEARMAN_lhrh_n6_net17_mrna_corrvals.csv'), avg.cort.striat.express.df)
write.csv(file=paste0(caret.dir, '/lhrh_n6_net17_rsfcmri_corrvals.csv'), striat2cort.fcmri)
reg.names <- colnames(avg.cort.striat.express.df)
for (reg in reg.names){
  write.vals   <- avg.cort.striat.express.df[[reg]]
  write.names  <- rownames(avg.cort.striat.express.df)
  write.idxs   <- 1:length(write.names)
  write.matrix <- cbind(write.idxs, write.vals, write.names)
  colnames(write.matrix) <- c('reg_idx','cor_vals','reg_name')
  write.csv(x=write.matrix, file=paste0(caret.dir, '/lhrh_mrna_net17_n6_striato_cort_', reg, '.csv'))
}
# -----------------------







# --------------------------------------------
# --------------------------------------------
# HUMAN STRIATAL REPLICATION - NIH GTEx DATA
# --------------------------------------------
# --------------------------------------------

# !!!
# !!! Please run the GTEx_preprocessing.R script first!
# !!!

# Read GTEx data
# ------------------
load(file=paste(base.dir, 'data/GTEx/nacc.contrast.Rdata', sep = '')) # NAcc vs Caudate+Putamen in GTEx data
load(file=paste(base.dir, 'data/GTEx/dge.counts.Rdata', sep = ''))
load(paste(base.dir, 'data/GTEx/ensembl2entrez.Rdata', sep = ''))
# ------------------


# Identify genes that are present in both the GTEx and AHBA datasets
# -------------
# entrez ids of genes wihtin the AHBA analysis
entrez.arr     <- all_data[[1]]$probes_collapse$entrez_id[order(all_data[[1]]$probes_collapse$gene_symbol)] # put the entrez ID in alphabetical order of HGNC symbols

# match the AHBA/GTEx sets based on entrez IDs. 
ref.entrezGene <- cbind(as.character(entrez.arr), as.character(sort(all_data[[1]]$probes_collapse$gene_symbol))) # columns containing entrez IDs and gene names
num.matches    <- length(which(rownames(nacc.contrast) %in% entrez.arr)) # ensembl2gene = GTEx genes
genes.in.both  <- rownames(nacc.contrast)[which(rownames(nacc.contrast) %in% entrez.arr)]
gtex.nacc.foldchange <- nacc.contrast[rownames(nacc.contrast) %in% genes.in.both,] # subset the GTEx data for matchs. 
gtex.nacc.foldchange <- gtex.nacc.foldchange[order(rownames(gtex.nacc.foldchange)),]
# -------------

# Convert the GTEx Entrez IDs to gene names (defined by AHBA, since HGNC symbols can vary depending on genome build)
# -------------
ref.entrezGene.tmp    <- ref.entrezGene[ref.entrezGene[,1] %in% rownames(gtex.nacc.foldchange),]
ref.entrezGene.sorted <- ref.entrezGene.tmp[order(ref.entrezGene.tmp[,1]),]
# check that the gene order was maintained...
check.matches <- length(which(ref.entrezGene.sorted[,1] == rownames(gtex.nacc.foldchange))) # all true if they are in the same order. 
if (check.matches != num.matches){
  print('GTEx and AHBA datasets did not cross-reference properly')
} else {
  print('all good')
}
rownames(gtex.nacc.foldchange) <- ref.entrezGene.sorted[,2]
# -------------

# create a list of GTEx significantly differentially expressed genes
pos.q01.idxs <- intersect(which(gtex.nacc.foldchange$adj.P.Val <= .01), which(gtex.nacc.foldchange$logFC > 0))
sig.nacc.DF  <- gtex.nacc.foldchange[pos.q01.idxs,]

# subset of matching genes between GTEx and AHBA
dge.counts.sub <- dge.counts
dge.counts.sub$counts <- dge.counts$counts[rownames(dge.counts$counts) %in% entrez.arr,]
dge.counts.sub$genes  <- dge.counts$genes$genes[dge.counts$genes$genes %in% entrez.arr]

lcpm <- cpm(dge.counts, log=TRUE)
lcpm.mean.norm <- t(apply(lcpm, 1, function(x) x-mean(x) ))

ref.entrezGene.sorted[ref.entrezGene.sorted[,2] == 'PVALB']
ref.entrezGene.sorted[ref.entrezGene.sorted[,2] == 'SSTR1']

# Plot the GTEx expression values
# ----------------------
regions2plot <- c('Nucleus_accumbens', 'Caudate', 'Putamen')
plot.gtex(lcpm.mean.norm, dge.counts.sub$samples$region, '5816', regions2plot, ylimit=3) # PVALB
plot.gtex(lcpm.mean.norm, dge.counts.sub$samples$region, '6751', regions2plot, ylimit=.5) # SSTR1


# correlate the fold change of limbic striatum AHBA to limbic striatum GTEx
# ----------------------
gtex.nacc.foldchange.sort <- gtex.nacc.foldchange[order(rownames(gtex.nacc.foldchange)),]
ahba.limbic.subset        <- striat.foldchanges.n6$stats$Limbic[rownames(striat.foldchanges.n6$stats$Limbic) %in% rownames(gtex.nacc.foldchange),]
length(rownames(ahba.limbic.subset) == rownames(gtex.nacc.foldchange))
cor.test(ahba.limbic.subset$logFC, gtex.nacc.foldchange.sort$logFC)


# FIGURE 6; Get significant AHBA genes, plot them against GTEx
# ----------------------
ahba.sig.limbic.genes <- striat.foldchanges.n6$q01_genes$Limbic
replicationPlot(ahba.array = ahba.limbic.subset$logFC, 
                ahba.sig.genes = ahba.sig.limbic.genes,
                rep.array = gtex.nacc.foldchange.sort$logFC,
                rep.fit = gtex.nacc.foldchange.sort, 
                xlabel = 'AHBA Limbic Log2 Fold Change',
                ylabel = 'GTEx NAcc (Limbic) Log Fold Change',
                xlimit = 5, 
                ylimit = 5)

write.csv(file=paste0(base.dir, '/output_files/GTEx_striatal_fold_change.csv'), x=gtex.nacc.foldchange.sort, quote=F)

                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
                        

# ---------------------------------------
# ---------------------------------------
# Brainspan Cortical Replication
# ---------------------------------------
# ---------------------------------------
brainspan_dir <- paste(base.dir, 'data/BRAINSPAN', sep = '')
gene.columns  <- read.csv(paste(brainspan_dir, '/columns_metadata.csv', sep = ''))
gene.rows     <- read.csv(paste(brainspan_dir, '/rows_metadata.csv', sep = ''))
expr.matrix   <- read.csv(paste(brainspan_dir, '/expression_matrix.csv', sep = ''), header = FALSE, row.names = 1)


# Identify overlapping genes between AHBA/Brainspan
# -----------------------
brainspan.genes.use <- which(gene.rows$entrez_id %in% all_data[[1]]$probes_collapse$entrez_id)
expr.matrix.use     <- expr.matrix[brainspan.genes.use,]
gene.rows.use       <- gene.rows[brainspan.genes.use,]
# -----------------------
# Dictionary with ensembl/entrez gene information
# -----------------------
human      <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- c("ensembl_gene_id", "entrezgene")
ensembl2entrez <- getBM( attributes, filters="entrezgene", values = all_data[[1]]$raw_probes$entrez_id, mart = human, bmHeader=FALSE)
# -----------------------

# Select adult donors
# -----------------------
adult.idxs   <- grep('18 yrs|19 yrs|21 yrs|23 yrs|30 yrs|36 yrs|37 yrs|40 yrs', gene.columns$age)
adult.expr   <- expr.matrix.use[,adult.idxs]
adult.cols   <- gene.columns[adult.idxs,] 
adult.donors <- unique(adult.cols$donor_name)
# -----------------------

# collapse brainspan data seperately for each individual subject
# -----------------------
rowID    <- rownames(adult.expr)
rowGroup <- gene.rows.use$entrez_id
all_collapsed_data    <- NULL
all_collapsed_columns <- NULL
donors                <- levels(droplevels(adult.cols$donor_name))
for ( donor in as.character(adult.donors) ) {
  print(donor)
  cur_expr          <- adult.expr[, adult.cols$donor_name == donor]
  log2_expr         <- log2(cur_expr + 1)
  cur_columns       <- adult.cols[adult.cols$donor_name == donor, ]
  collapsed_matrix  <- collapseRows(log2_expr, rowGroup, rowID, 
                                    method="MaxMean", connectivityBasedCollapsing=FALSE,
                                    connectivityPower=1, selectFewestMissing=TRUE, thresholdCombine=NA)
  collapsed_data        <- collapsed_matrix$datETcollapsed
  all_collapsed_data    <- cbind(all_collapsed_data, collapsed_data)
  all_collapsed_columns <- rbind(all_collapsed_columns, cur_columns)
}
# save for quicker reading later'
save(file = paste(brainspan_dir, '/adult_columns_metadata.Rdata', sep = ''), x = all_collapsed_columns)
save(file = paste(brainspan_dir, '/adult_expression_matrix.Rdata', sep = ''), x = all_collapsed_data)
load(file = paste(brainspan_dir, '/adult_columns_metadata.Rdata', sep = ''))
load(file = paste(brainspan_dir, '/adult_expression_matrix.Rdata', sep = ''))


# Cortical regions for in the analysis
# -------------------
use_regs <- c('OFC','M1C','V1C','A1C','ITC','VFC','MFC','STC','DFC','S1C','IPC')
use_cols <- all_collapsed_columns[all_collapsed_columns$structure_acronym %in% use_regs,]
# -----------------------

# Map the Brainspan entrez IDs onto AHBA HGNC gene symbols
# -------------------
use_idxs   <- which(all_collapsed_columns$structure_acronym %in% use_regs)
use_cols   <- all_collapsed_columns[use_idxs,]
bspan.eset <- all_collapsed_data[,use_idxs]
colnames(bspan.eset) <- use_cols$structure_acronym
bspan.eset <- bspan.eset[which(rowSums(bspan.eset) != 0),] # git rid of genes with no exprssion 

ahba.bspan.subset <- all_data[[1]]$probes_collapse[which(all_data[[1]]$probes_collapse$entrez_id %in% rownames(bspan.eset)),]

# sort ahba by entrez id 
ahba.bspan.sort <- ahba.bspan.subset[order(as.character(ahba.bspan.subset$entrez_id)),]
length(which(ahba.bspan.sort$entrez_id == rownames(bspan.eset)))
dim(ahba.bspan.subset)
dim(bspan.eset)
rownames(bspan.eset) <- ahba.bspan.sort$gene_symbol

# -----------------------
# Significant network-associated AHBA genes 
# -------------------
ahba.sig.limbic.genes  <- intersect(striat.foldchanges.n6$q01_genes$Limbic, cort.foldchanges.n6$q01_genes$Limbic)
ahba.sig.sommot.genes  <- intersect(striat.foldchanges.n6$q01_genes$SomMot, cort.foldchanges.n6$q01_genes$SomMot)
# -----------------------

# Plot expression across the brainspan adult data
# -------------------
cort.regions2plot <- c('OFC','MFC','ITC','VFC','DFC','IPC','STC','M1C','A1C','S1C','V1C')
plot.brainspan(bspan.eset, use_cols, 'PVALB', cort.regions2plot, ylimit=2.2)
plot.brainspan(bspan.eset, use_cols, 'SSTR1', cort.regions2plot, ylimit=1.5)
# -------------------

# Create the design/contrast matrix
# -------------------
use_cols         <- droplevels(use_cols)
use_cols$age_num <- as.numeric(as.character(gsub(' yrs', '', use_cols$age)))
design           <- model.matrix(~0 + structure_acronym + age_num + gender, data=use_cols)
colnames(design) <- gsub('structure_acronym', '', colnames(design))
donor_array      <- use_cols$donor_name
use_cols         <- droplevels(use_cols)

# -------------------
# SOMATOMOTOR vs REST OF CORTEX
# -------------------
cur_contrast <- "-.125*OFC+0.333333*M1C-0.125*V1C+0.333333*A1C-.125*ITC-0.125*VFC-.125*MFC-0.125*STC-0.125*DFC+0.333333*S1C-0.125*IPC"
cm           <- makeContrasts(cur_contrast, levels=design)
corfit       <- duplicateCorrelation(bspan.eset, design=design, block=donor_array)
fit_brsp     <- lmFit(bspan.eset, design, block=donor_array, correlation=corfit$consensus.correlation)
fit.brsp.sommot <- eBayes( contrasts.fit( fit_brsp, cm ), trend=TRUE)
sommot.dex   <- topTable(fit.brsp.sommot, number = Inf)
sommot.stats <- sommot.dex[order(rownames(sommot.dex)),]
gene.rows.ordered <- gene.rows.use[order(gene.rows.use$gene_symbol),]
plot.ahba.sommot  <- cort.foldchanges.n6$stats$SomMot

ahba.sommot.subset <- plot.ahba.sommot[rownames(plot.ahba.sommot) %in% rownames(sommot.stats),]
cor.test(ahba.sommot.subset$logFC, sommot.stats$logFC)
write.csv(x=sommot.stats[order(rownames(sommot.stats)),], file=paste0(base.dir, 'output_files/Brainspan_sommot_cortical_foldchange.csv'))

# ----------------------
# Significt somato/motor Brainspan genes
# ----------------------
sig_idxs  <- which(sommot.stats$adj.P.Val <= 0.01)
pos_idxs  <- which(sommot.stats$logFC > 0)
use.idxs  <- intersect(sig_idxs, pos_idxs)
sommot.idxs     <- which(all_data[[1]]$probes_collapse$gene_symbol %in% cort.foldchanges.n6$q01_genes$SomMot)
sommot.sig.ahba <- all_data[[1]]$probes_collapse$entrez_id[sommot.idxs]
# ----------------------

glength   <- dim(cort.foldchanges.n6$stats$Default)[1]
ahba.idxs <- 1:glength
replicationPlot(ahba.array = ahba.sommot.subset$logFC, 
                ahba.sig.genes = cort.foldchanges.n6$q01_genes$SomMot,
                rep.array = sommot.stats$logFC,
                rep.fit = sommot.stats, 
                xlabel = 'AHBA SomMot Cortex Log2 Fold Change (n=6)',
                ylabel = 'Brainspan SomMot Log2 Fold Change (n=8)',
                xlimit = 2, 
                ylimit = 2)


# ----------------------
# Limbic vs REST OF CORTEX
# ----------------------
cur_contrast    <- "0.333333*OFC-0.125*M1C-0.125*V1C-0.125*A1C+0.333333*ITC-0.125*VFC+0.333333*MFC-0.125*STC-0.125*DFC-0.125*S1C-0.125*IPC"
cm              <- makeContrasts(cur_contrast, levels = design)
corfit          <- duplicateCorrelation(bspan.eset, design = design,  block = donor_array)
fit_brsp        <- lmFit(bspan.eset, design, block = donor_array, correlation=corfit$consensus.correlation)
#fit.brsp.limbic <- eBayes( contrasts.fit( fit_brsp, cm ))
fit.brsp.limbic <- eBayes( contrasts.fit( fit_brsp, cm ), trend=TRUE)


# ----------------------
limbic.dex   <- topTable(fit.brsp.limbic, number=Inf)
limbic.stats <- limbic.dex[order(rownames(limbic.dex)),]
# ----------------------

gene.rows.ordered <- gene.rows.use[order(gene.rows.use$gene_symbol),]
plot.ahba.limbic  <- cort.foldchanges.n6$stats$Limbic
#
ahba.limbic.subset <- plot.ahba.limbic[rownames(plot.ahba.limbic) %in% rownames(limbic.stats),]
cor.test(ahba.limbic.subset$logFC, limbic.stats$logFC)
write.csv(x=limbic.stats[order(rownames(limbic.stats)),], file=paste0(base.dir, 'output_files/Brainspan_limbic_cortical_foldchange.csv'))

# ----------------------
# Significt somato/motor Brainspan genes
# ----------------------
sig_idxs  <- which(limbic.stats$adj.P.Val <= 0.01)
pos_idxs  <- which(limbic.stats$logFC > 0)
use.idxs  <- intersect(sig_idxs, pos_idxs)
limbic.idxs     <- which(all_data[[1]]$probes_collapse$gene_symbol %in% cort.foldchanges.n6$q01_genes$Limbic)
limbic.sig.ahba <- all_data[[1]]$probes_collapse$entrez_id[limbic.idxs]
# ----------------------

ahba.idxs <- 1:glength
glength   <- dim(cort.foldchanges.n6$stats$Default)[1]
replicationPlot(ahba.array = ahba.limbic.subset$logFC, 
                ahba.sig.genes = cort.foldchanges.n6$q01_genes$Limbic,
                rep.array = limbic.stats$logFC,
                rep.fit = limbic.stats, 
                xlabel = 'AHBA Limbic Cortex Log2 Fold Change (n=6)',
                ylabel = 'Brainspan Limbic Log2 Fold Change (n=8)',
                xlimit = 2.5, 
                ylimit = 2.5)


# Genes significantly differentially expressed in Brainspan Som/Mot Cortex
# ------------------
pos.idxs <- which(sommot.stats$logFC > 0)
sig.idxs <- which(sommot.stats$adj.P.Val <= 0.01)
use.idxs <- intersect(pos.idxs, sig.idxs)
significant.cort.sommot <- rownames(sommot.stats)[use.idxs]

# Genes significantly differentially expressed in Brainspan Limbic Cortex
# ------------------
pos.idxs <- which(limbic.stats$logFC > 0)
sig.idxs <- which(limbic.stats$adj.P.Val <= 0.01)
use.idxs <- intersect(pos.idxs, sig.idxs)
significant.cort.limbic <- rownames(limbic.stats)[use.idxs]
# ----------------------
# Genes significantly differentially exprssion in both Brainspan & AHBA limbic network
# ------------------ 
pos.idxs <- which(gtex.nacc.foldchange.sort$logFC > 0)
sig.idxs <- which(gtex.nacc.foldchange.sort$adj.P.Val <= 0.01)
use.idxs <- intersect(pos.idxs, sig.idxs)
significant.striat.limbic <- rownames(gtex.nacc.foldchange.sort)[use.idxs] # replication GTEx striatal genes
replication.limbic.genes  <- intersect(significant.cort.limbic, significant.striat.limbic)
ahba.limbic.genes    <- intersect(striat.foldchanges.n6$q01_genes$Limbic, cort.foldchanges.n6$q01_genes$Limbic)
sig.in.both.datasets <- intersect(ahba.limbic.genes, replication.limbic.genes)

# ----------------------
# Hypergeometric tests
# Determine if there is significant amounts of overlap between GTEx LIMBIC STR to BSPAN LIMBIC CORTEX
# ------------------ 
all.bspan.genes    <- rownames(limbic.stats) # Number of Brainspan genes
num.overlap        <- length(replication.limbic.genes) # number of GTEx+Brainspan genes
num.genes          <- length(which(rownames(gtex.nacc.foldchange.sort) %in% all.bspan.genes))
white.genes        <- length(significant.striat.limbic)
black.genes        <- num.genes-length(significant.striat.limbic)
genes.drawn        <- length(significant.cort.limbic)
print(num.overlap)
phyper(num.overlap-1, white.genes, black.genes, genes.drawn, lower.tail = F)

# ----------------------
# Determine if there is significant amounts of overlap between limbic associated genes in replication data
# ------------------ 
replication.genes  <- rownames(gtex.nacc.foldchange)[which(rownames(gtex.nacc.foldchange) %in% all.bspan.genes)] # number of genes that overlap in GTEx/BRAINSPAN
num.genes          <- length(which(replication.genes %in% all_data[[1]]$probes_collapse$gene_symbol)) # genes common across all datasets
sig.in.both.datasets  <- intersect(ahba.limbic.genes, replication.limbic.genes) # significant overlapping limbic genes
num.overlap        <- length(sig.in.both.datasets)
white.genes        <- length(replication.limbic.genes)
black.genes        <- num.genes-length(replication.limbic.genes)
genes.drawn        <- length(ahba.limbic.genes)
print(num.overlap)
phyper(num.overlap-1, white.genes, black.genes, genes.drawn, lower.tail = F)
# ----------------------



                          
                          
# -----------------------------
# -----------------------------
# Compare limbic cortico-striatal genes to with previously identified gene sets (e.g. Richiardi 2015; Wang 2015)
# -----------------------------
# -----------------------------
ref.gene.sets <- read.csv(paste(base.dir, 'reference_files/Richiardi_Wang_Gene_Sets.csv', sep = ''))

# ----------------------
# Overlap and hypergeometric test for Richiardi Gene Set
# ----------------------
num.richiardi <- which(ref.gene.sets$Richiardi_Genes %in% ahba.limbic.genes)
num.overlap   <- length(num.richiardi)
num.genes     <- length(cort.foldchanges.n6$fit_df$Limbic$coefficients)
white.genes   <- length(ahba.limbic.genes)
black.genes   <- num.genes-length(ahba.limbic.genes)
genes.drawn   <- length(ref.gene.sets$Richiardi_Genes)
phyper(num.overlap-1, white.genes, black.genes, genes.drawn, lower.tail = F)
# -------------
# Overlap and hypergeometric test for Wang Gene Set
# -------------
num.wang     <- which(ref.gene.sets$Wang_Genes %in% ahba.limbic.genes)
num.overlap  <- length(num.wang)
num.genes    <- length(cort.foldchanges.n6$fit_df$Limbic$coefficients)
white.genes  <- length(ahba.limbic.genes)
black.genes  <- num.genes-length(ahba.limbic.genes)
genes.drawn  <- length(ref.gene.sets$Wang_Genes)
phyper(num.overlap-1, white.genes, black.genes, genes.drawn, lower.tail = F)
# -------------


                          
                          

# -----------------------------
# -----------------------------
# Identify cell type enrichment of the limbic-network genes according to the Barres data
# -----------------------------
# -----------------------------
# Create BiomaRt data dictionaries to convert gene IDs between human/mouse
# Map Human ensembl to mouse ensembl
# -------------------
human      <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_associated_gene_name")
orth.mouse <- getBM( attributes, filters="with_mmusculus_homolog", values = TRUE, mart = human, bmHeader=FALSE)
# ----------------------

# ----------------------
# Human ensembl to symbol and entrez ID
# ----------------------
entrez2ensemble <- getBM(attributes = c('entrezgene', 'hgnc_symbol', 'entrezgene_trans_name', 'ensembl_gene_id'), values = TRUE, mart = human, bmHeader=FALSE)
entrez.dict     <- entrez2ensemble[!is.na(entrez2ensemble$entrezgene),]

# ----------------------
# Mouse Affymetrix IDs to ensemble IDs
# -------------------
mouse          <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
attributes     <- c("affy_mouse430_2","ensembl_gene_id")
mouse.affy2ens <- getBM(attributes = attributes, values = TRUE, mart = mouse, bmHeader=FALSE)
use.mouse.affy2ens <- mouse.affy2ens[mouse.affy2ens$affy_mouse430_2 != '',]

# ----------------------
# Ensembl/Entrez ID of limbic cortico-striatal genes
# ----------------------
entrez.limbic    <- all_data[[1]]$probes_collapse$entrez_id[which(all_data[[1]]$probes_collapse$gene_symbol %in% limbic.genes)]
ensembl.limbic   <- entrez.dict$ensembl_gene_id[which(entrez.dict$entrezgene %in% entrez.limbic)]
mouse.gene.names <- orth.mouse$mmusculus_homolog_associated_gene_name[orth.mouse$ensembl_gene_id %in% ensembl.limbic]
mouse.gene.names <- unique(mouse.gene.names)
# ----------------------

# Create a list of human genes with identifable mouse orthologs
# ----------------------
human2mouse <- NULL
for (m.gene in mouse.gene.names){
  hum.ensembl <- orth.mouse$ensembl_gene_id[orth.mouse$mmusculus_homolog_associated_gene_name == m.gene]
  hum.genes   <- NULL
  for (ens in hum.ensembl){
    ent <- entrez.dict$entrezgene[entrez.dict$ensembl_gene_id == ens]
    #print(ent)
    hum.gene <- as.character(all_data[[1]]$probes_collapse$gene_symbol[which(all_data[[1]]$probes_collapse$entrez_id == ent)])
    hum.genes <- c(hum.genes, hum.gene)
  }
  hum.genes <- unique(hum.genes)
  hum.genes <- hum.genes[hum.genes %in% limbic.genes]
  human2mouse <- c(human2mouse, paste(hum.genes, collapse = '|'))
}
head(cbind(mouse.gene.names,human2mouse)) # inspect the gene mapping
# ----------------------


# Examine whether limbic cortico-striatal genes are expressed most in one cell type 
# ------------------ 
cell.type.dat    <- read_excel(paste(base.dir, 'data/BARRES_CELLTYPE/barreslab_rnaseq.xlsx', sep = ''))
enrichment.array <- NULL
fc.array         <- NULL
fc.array.strict  <- NULL
gene.array       <- NULL
human.genes      <- NULL
ct <- 1
for (gene in mouse.gene.names){
  human.gene <- human2mouse[ct]
  ct <- ct + 1
  gene.idx <- which(toupper(cell.type.dat$`Gene symbol`) %in% toupper(gene)) # index of this gene in the Barres data
  gene.row <- cell.type.dat[gene.idx,]
  ncols    <- dim(cell.type.dat)[2]
  gene.rows.numbers <- gene.row[1,3:ncols]
  max.col  <- which(gene.rows.numbers == max(gene.rows.numbers))
  
  # Identify which cell type has highest expression
  if (length(max.col) != 0){
    idxs <- 1:7
    comparison.idxs <- idxs[!idxs %in% max.col]
    max.cell        <- names(gene.rows.numbers)[max.col] # cell type with the highest expression
    if (length(max.cell) == 1){
      # compare the expression of the max cell type to the expression of the next highest cell type
      foldchange.strict <- log2(as.numeric(gene.rows.numbers[max.col])) - log2(max(gene.rows.numbers[comparison.idxs]))
      fc.array.strict   <- c(fc.array.strict, foldchange.strict) 
      enrichment.array  <- c(enrichment.array, max.cell)
      gene.array        <- c(gene.array, gene)
      human.genes <- c(human.genes, human.gene)
    }
  } else { # not all genes are present in the Barres data
    print(gene)
  }
}
# Construct a data frame with the enrichment data output
# ------------------ 
enrich.out                 <- as.data.frame(cbind(gene.array, enrichment.array, fc.array, fc.array.strict, human.genes))
enrich.out$fc.array.strict <- as.numeric(as.character(enrich.out$fc.array.strict))
enrich.out$fc.array <- as.numeric(as.character(enrich.out$fc.array))


# Number of genes that show 1.5 enrichment in any category
# ------------------ 
enrich.out.sig <- enrich.out[which(enrich.out$fc.array.strict > 1.5),]
cell.table     <- rev(sort(table(enrich.out[which(enrich.out$fc.array.strict > 1.5),]$enrichment.array)))
# ----------------------
# Format and write the output
# ------------------ 
enrich.genes.out <- vector()
names <- NULL
nums  <- NULL
for (name in names(cell.table)){
  genes.tmp        <- as.character(enrich.out.sig[enrich.out.sig$enrichment.array == name,]$gene.array)
  enrich.genes.out <- rbind.fill.matrix(enrich.genes.out, t(genes.tmp))
  genes.tmp        <- as.character(enrich.out.sig[enrich.out.sig$enrichment.array == name,]$human.genes)
  enrich.genes.out <- rbind.fill.matrix(enrich.genes.out, t(genes.tmp))
  names <- c(names, rep(name,2))
  nums  <- c(nums, rep(cell.table[[name]], 2))
}
enrich.genes.out[is.na(enrich.genes.out)] <- ''
t.write.genes <- t(enrich.genes.out)
colnames(t.write.genes) <- names
write.me <- rbind(names, nums, t.write.genes)
write.csv(file = paste(base.dir, 'output_files/barres_cell_enrichment_counts.csv', sep = ''), x=write.me)
# ----------------------




                          
                          
                          
                          
                          
# -------------------------------------------
# -------------------------------------------
# Striatal NIH Blueprint Non-human primate data 
# -------------------------------------------
# -------------------------------------------
blueprint.dir     <- paste(base.dir, 'data/BLUEPRINT/', sep = '')
blueprint.columns <- read.csv(file = paste(blueprint.dir, 'columns_metadata.csv', sep = ''))
blueprint.rows    <- read.csv(file = paste(blueprint.dir, 'rows_metadata.csv', sep = ''))
expression_data   <- fread(paste(blueprint.dir, 'expression_matrix.csv', sep=''), header=FALSE, drop=1)

# save(file = paste(blueprint.dir, 'expression_matrix.Rdata', sep = ''), x=expression_data)
load(file = paste(blueprint.dir, 'expression_matrix.Rdata', sep = ''))
expression_data <- as.matrix(expression_data)
colnames(expression_data) <- blueprint.columns$structure_acronym
rownames(expression_data) <- blueprint.rows$probe_name
# ------------------------

# ----------------------
# Only analyze high probability human/macaque gene homologs identified by Bakken et al., 2016
# ------------------------
bakken.gene.mapping <- read.csv(paste(base.dir, 'data/BAKKEN/bakken_probe_mapping.csv', sep = ''))
valid.human.entrez  <- which(!is.na(bakken.gene.mapping$human_entrezid)) # bakken probes with human entrez
keep.for.analysis   <- which(bakken.gene.mapping$keep_for_analysis) # probes kept for analyses by Bakken et al., 
use.these.idxs      <- intersect(valid.human.entrez, keep.for.analysis) 
probes.to.use       <- bakken.gene.mapping$probeid[use.these.idxs]
probes.to.use       <- droplevels(probes.to.use)
# ----------------------
# Subset the NIH Blueprint data according to the high-confidence probes
# ---------------------
blueprint.expr.filt <- expression_data[rownames(expression_data) %in% probes.to.use,]
blueprint.rows.filt <- blueprint.rows[blueprint.rows$probe_name %in% probes.to.use,]
# ----------------------

# Convert the probe names into human gene symbols/entrez ids
# ------------------------
# subset bakken and blueprint datasets based on matching probeids
bakken.subset       <- bakken.gene.mapping[which(bakken.gene.mapping$probeid %in% blueprint.rows.filt$probe_name),]
bakken.subset.sort  <- bakken.subset[order(bakken.subset$probeid),]
blueprint.rows.sort <- blueprint.rows.filt[order(blueprint.rows.filt$probe_name),]
length(which(blueprint.rows.sort$probe_name == bakken.subset.sort$probeid))

blueprint.rows.filt.wHomologs <- cbind(cbind(blueprint.rows.sort, as.character(bakken.subset.sort$human_genesymbol)), as.character(bakken.subset.sort$human_entrezid))
colnames(blueprint.rows.filt.wHomologs) <- c(colnames(blueprint.rows.filt), 'humanGene', 'humanEntrez')


# note: there is not a 1-to-1 mapping between each probe and each human gene.
# Remove small subset of probes that have non-unique mapping to the human genome
gene.counts     <- table(blueprint.rows.filt.wHomologs$humanGene)
nonunique_genes <- names(gene.counts[gene.counts > 1])
blueprint.rows.unique.wHomologs <- blueprint.rows.filt.wHomologs[!blueprint.rows.filt.wHomologs$humanGene %in% nonunique_genes,]
# ----------------------

# Subset the expression data
# ------------------------
blueprint.eset.use <- blueprint.expr.filt[rownames(blueprint.expr.filt) %in% blueprint.rows.unique.wHomologs$probe_name ,]
blueprint.eset.use <- blueprint.eset.use[order(rownames(blueprint.eset.use)),] # order expression data
blueprint.rows.use <- blueprint.rows.filt.wHomologs[blueprint.rows.filt.wHomologs$probe_name %in% blueprint.rows.unique.wHomologs$probe_name, ]
blueprint.rows.use <- blueprint.rows.use[order(blueprint.rows.use$probe_name),] # order row data
length(which(blueprint.rows.use$probe_name == rownames(blueprint.eset.use)))
# ----------------------
# Select the striatal column indices
# ------------------------
striatal.reg       <- c('Ca', 'NAC', 'Pu')
striatal.idxs      <- sort(which(blueprint.columns$structure_acronym %in% striatal.reg))
# ------------------------
# ----------------------
# Subset Striatal Data
# ------------------------
blueprint.striat.expr <- blueprint.eset.use[,striatal.idxs]
blueprint.striat.cols <- blueprint.columns[striatal.idxs,]
# ------------------------


# Subset adolescent/adult macaque brains
# ------------------------
mature.brains          <- c(grep('12 mo', blueprint.striat.cols$age), grep('48 mo', blueprint.striat.cols$age))
blueprint.mature.cols  <- blueprint.striat.cols[mature.brains,]
blueprint.mature.expr  <- blueprint.striat.expr[, mature.brains]
blueprint.mature.cols  <- blueprint.mature.cols[which(!is.na(blueprint.mature.expr[1,])),]
blueprint.mature.expr  <- blueprint.mature.expr[,which(!is.na(blueprint.mature.expr[1,]))]
# ------------------------
# Formate donor and structure information
# ------------------------
donors        <- unique(blueprint.mature.cols$donor_name)
donor.arr     <- blueprint.mature.cols$donor_name
reg.arr       <- blueprint.mature.cols$structure_acronym
regions       <- c('NAC','Ca','Pu')
blueprint.mature.cols$age_num <- as.numeric(as.character(gsub(' mo', '', blueprint.mature.cols$age)))
# ------------------------

rownames(blueprint.mature.expr) <- blueprint.rows.use$humanGene

blueprint.mature.cols <- droplevels(blueprint.mature.cols)
design           <- model.matrix(~0 + structure_acronym + age_num, data=blueprint.mature.cols)
colnames(design) <- gsub('structure_acronym', '', colnames(design))
corfit           <- duplicateCorrelation(blueprint.mature.expr, design, block=blueprint.mature.cols$donor_name)
cur.contrast     <- '1*NAC-.5*Ca-.5*Pu'
cmtx             <- makeContrasts(cur.contrast, levels = design) # Make a contrast matrix
tmplm            <- lmFit(blueprint.mature.expr, design, block=donor.arr, correlation=corfit$consensus.correlation ) # Fit the model to the data
fit.blueprint.nacc  <- eBayes(contrasts.fit(tmplm, cmtx))
blueprint.nacc.table <- topTable(fit.blueprint.nacc, num = Inf)
blueprint.nacc.table.sort <- blueprint.nacc.table[order(rownames(blueprint.nacc.table)),]
  
regions2plot <- c('NAC', 'Ca', 'Pu')
plot.blueprint(blueprint.mature.expr, colnames(blueprint.mature.expr), 'SSTR1', regions2plot, ylimit=2)


# Cross-reference to the human gtex and ahba striatal data
# ------------------------
blueprint.rows.sort  <- blueprint.rows.use[order(blueprint.rows.use$humanGene),]
rownames(blueprint.nacc.table.sort) <- as.character(blueprint.rows.sort$humanEntrez)

ahba.probe.info      <- all_data[[1]]$probes_collapse
ahba.probe.info.sort <- ahba.probe.info[order(ahba.probe.info$gene_symbol),] # sort by gene symbol
length(which(ahba.probe.info.sort$gene_symbol == rownames(striat.foldchanges.n6$stats$Limbic))) # check order

# probes in both AHBA and blueprint data tables
ahba.entrez.matches        <- droplevels(blueprint.rows.sort$humanEntrez[which(rownames(blueprint.nacc.table.sort) %in% ahba.probe.info.sort$entrez_id)])
blueprint.nacc.ahba.subset <- blueprint.nacc.table.sort[rownames(blueprint.nacc.table.sort) %in% ahba.entrez.matches,]
blueprint.nacc.ahba.subset <- blueprint.nacc.ahba.subset[order(rownames(blueprint.nacc.ahba.subset)),]

limbic.sortByEntrez <- striat.foldchanges.n6$stats$Limbic[order(as.character(ahba.probe.info.sort$entrez_id)),]
rownames(limbic.sortByEntrez) <- sort(as.character(ahba.probe.info.sort$entrez_id))
ahba.limbic.subsetByBlueprint <- limbic.sortByEntrez[rownames(limbic.sortByEntrez) %in% rownames(blueprint.nacc.ahba.subset),]

# BLUEPRINT - AHBA limbic striatum fold change correlation 
# ---------------
cor.test(ahba.limbic.subsetByBlueprint$logFC, blueprint.nacc.ahba.subset$logFC)

# BLUEPRINT - GTEx limbic striatum fold change correlation 
# ---------------
gtex.nacc.subsetByBlueprint <- gtex.nacc.foldchange[gtex.nacc.foldchange$genes %in% rownames(blueprint.nacc.table.sort),]
blueprint.nacc.subsetByGTEx <- blueprint.nacc.table.sort[rownames(blueprint.nacc.table.sort) %in% gtex.nacc.foldchange$genes,]
blueprint.nacc.subsetByGTEx <- blueprint.nacc.subsetByGTEx[order(rownames(blueprint.nacc.subsetByGTEx)),]
cor.test(gtex.nacc.subsetByBlueprint$logFC, blueprint.nacc.subsetByGTEx$logFC)

# ------------------------
# Plot 
# ------------------------
ahba.nacc.entrez <- rownames(ahba.limbic.subsetByBlueprint)[intersect(which(ahba.limbic.subsetByBlueprint$logFC > 0), which(ahba.limbic.subsetByBlueprint$adj.P.Val <= .01))]
write.csv(x=blueprint.nacc.ahba.subset, file=paste0(base.dir, '/output_files/Blueprint_NAcc_DEX.csv'))
monk.replicationPlot(ahba.dat = ahba.limbic.subsetByBlueprint, 
                     ahba.sig.genes = ahba.nacc.entrez,
                     rep.dat = blueprint.nacc.ahba.subset,
                     xlabel = 'AHBA Limbic striatal Log2 Fold Change (n=6)',
                     ylabel = 'Blueprint Macaque NACC',
                     xlimit = 5, 
                     ylimit = 5)


                          
                          
                          
                          
                          

# --------------------------------------------------
# --------------------------------------------------
# Bernard Cortical Data
# --------------------------------------------------
# --------------------------------------------------
# RMA and batch normalized log10 expression data is available on NIH GSE Data website
bernard.dir  <- paste(base.dir, 'data/BERNARD', sep = '')
gse2         <- getGEO(filename=paste0(bernard.dir, '/GSE31613_series_matrix.txt'), 
                       GSEMatrix=TRUE)
# save(file=paste(bernard.dir, '/GSE31613_series_matrix.Rdata', sep = ''), x=gse2)          
load(file=paste(bernard.dir, '/GSE31613_series_matrix.Rdata', sep = ''))        
bernard.gse <- gse2

# ----------------------
# Extract the expression values, make sure they are in log2 form
# -----------------
bernard.expr.matrix <- bernard.gse@assayData$exprs
bernard.expr.log2   <- log2(10^bernard.expr.matrix) # convert from log10 to log2 
colnames(bernard.expr.log2) <- bernard.gse@phenoData@data$source_name_ch1 # set region names to column names
# -----------------
# Define the regions corresponding to limbic cortex
# -----------------
limbic.cortex    <- c('Anterior cingulate gyrus','Orbitofrontal cortex')
nonlimbic.cortex <- c('Primary auditory ctx','Dorsolateral Prefrontal Cortex','Primary somatosensory cortex',
                      'Primary motor cortex', 'Temporal area', 'middle temporal area', 'Middle temporal area','V1', 'V2')
bernard.regions  <- c(limbic.cortex, nonlimbic.cortex)
# -----------------

# Only analyze high probability human/macaque gene homologs identified by Bakken et al., 2016
# ------------------------
bakken.gene.mapping <- read.csv(paste(base.dir, 'data/BAKKEN/bakken_probe_mapping.csv', sep = ''))
valid.human.entrez  <- which(!is.na(bakken.gene.mapping$human_entrezid)) # bakken probes with human entrez
keep.for.analysis   <- which(bakken.gene.mapping$keep_for_analysis) # probes kept for analyses by Bakken et al., 
use.these.idxs      <- intersect(valid.human.entrez, keep.for.analysis) 
probes.to.use       <- bakken.gene.mapping$probeid[use.these.idxs]
# -----------------
# determine the indices of samples corresponding to cortical regions that we want to analyze
# ------------------------
use.idxs.tmp         <- lapply(bernard.regions, grep, x=colnames(bernard.expr.log2))
col.idxs             <- sort(unlist(use.idxs.tmp))
# -----------------

# Subset log2 expression data based on valid Bakken probes and cortical regions
# ------------------------
bernard.expr.subset <- bernard.expr.log2[which(rownames(bernard.expr.log2) %in% probes.to.use), col.idxs]
bernard.expr.sort   <- bernard.expr.subset[order(rownames(bernard.expr.subset)), ]

# Create a reference data frame with information about human/monkey entrez ids and gene names
# ------------------------
bakken.subset      <- bakken.gene.mapping[bakken.gene.mapping$probeid %in% rownames(bernard.expr.subset),]
bakken.subset.sort <- bakken.subset[order(bakken.subset$probeid),]
length(which(rownames(bernard.expr.sort) == bakken.subset.sort$probeid))


# Extract layer and subject information from the column names
# -------------------------
bernard.col.names <- colnames(bernard.expr.subset)
layer.array <- NULL
subj.arr    <- NULL
sex.arr     <- NULL
# extract layer and subject info from column names
for (col.name in bernard.col.names){
  layer       <- strsplit(strsplit(col.name, ',')[[1]][2], ';')[[1]][1]
  layer.name  <- gsub(' ','',layer)
  layer.array <- c(layer.array, layer.name)
  subj        <- strsplit(strsplit(col.name, ',')[[1]][2], ';')[[1]][2]
  subj.name   <- gsub(' ','',subj)
  subj.arr    <- c(subj.arr, subj.name)
}
sex.arr <- matrix(NA, length(subj.arr), 1)
sex.arr[grep('Male', subj.arr)] <- 1
sex.arr[grep('Female', subj.arr)] <- 0
# -------------------------
donors        <- unique(subj.arr)
donors        <- donors[!is.na(donors)]
donors <- gsub('Male|Female', '', donors)
averaged.data <- NULL
cols.out      <- NULL
donor.arr     <- NULL
reg.arr       <- NULL
sex.arr.collapsed <- NULL
# For each region
for (reg in bernard.regions){
  reg.idxs <- grep(reg, bernard.col.names)
  # for each donor
  for ( donor in donors ){
    don.idxs  <- grep(donor, bernard.col.names)
    cur.idxs  <- intersect(don.idxs, reg.idxs)
    if (length(cur.idxs) != 0){
      donor.arr <- c(donor.arr, donor)
      reg.arr   <- c(reg.arr, reg)
      sex.arr.collapsed <- c(sex.arr.collapsed, sex.arr[cur.idxs][1])
      # average data for this region/donor if more than one sample exists
      cur.don.reg   <- rowMeans(bernard.expr.subset[,cur.idxs])
      averaged.data <- cbind(averaged.data, cur.don.reg)
      cols.out      <- c(cols.out, paste(donor, reg, sep = ' '))
    }
  }
}
colnames(averaged.data) <- cols.out
averaged.data.use       <- averaged.data 


# get rid of genes that have non-unique correspondence between human/primates
# -------------------------
gene.counts     <- table(as.character(bakken.subset.sort$human_genesymbol))
nonunique_genes <- names(gene.counts[gene.counts > 1])
bernard.wHomologs.unique <- bakken.subset.sort[!bakken.subset.sort$human_genesymbol %in% nonunique_genes,]
averaged.data.unique     <- averaged.data.use[which(rownames(averaged.data.use) %in% bernard.wHomologs.unique$probeid),]
bernard.expr.unique      <- bernard.expr.subset[which(rownames(bernard.expr.subset) %in% bernard.wHomologs.unique$probeid), ]
colnames(averaged.data.unique) <- gsub(' ', '_', gsub('fac', '', colnames(averaged.data.unique)))

# Differential expresssino for limbic cortex vs. 8 other cortical regions
# -------------------------
colnames(averaged.data.unique) <- reg.arr
rownames(averaged.data.unique) <- as.character(bernard.wHomologs.unique$human_entrezid)
fac              <- as.factor(reg.arr)
sex.cov <- as.factor(sex.arr.collapsed)
design           <- model.matrix(~fac + sex.cov + 0)
colnames(design) <- gsub(' ', '_', gsub('fac', '', colnames(design)))
corfit           <- duplicateCorrelation(averaged.data.unique, design, block=donor.arr)
cur_contrast <- paste('.5*Orbitofrontal_cortex+.5*Anterior_cingulate_gyrus-0.125*Middle_temporal_area-',
                      '0.125*Dorsolateral_Prefrontal_Cortex-0.125*Primary_auditory_ctx-',
                      '0.125*Primary_somatosensory_cortex-0.125*V1-0.125*V2-0.125*Temporal_area-0.125*Primary_motor_cortex', sep = '')
cmtx         <- makeContrasts(cur_contrast, levels=design ) # Make a contrast matrix
tmplm        <- lmFit( averaged.data.unique, design, block=donor.arr, correlation=corfit$consensus.correlation ) # Fit the model to the data
bernard.limbic.fit.cort <- eBayes( contrasts.fit( tmplm, cmtx ) )
bernard.limbic.table    <- topTable(bernard.limbic.fit.cort, num=Inf)


# compare blueprint to AHBA
#rownames(bernard.limbic.table) <- rownames(bernard.limbic.table))
bernard.limbic.sort <- bernard.limbic.table[order(rownames(bernard.limbic.table)),]
write.csv(file=paste0(base.dir, '/output_files/bernard_cort_limbic.csv'), x=bernard.limbic.sort)


ahba.probes      <- all_data[[1]]$probes_collapse
ahba.probes.sort <- ahba.probes[order(ahba.probes$gene_symbol),]
cort.ahba.limbic <- cort.foldchanges.n6$stats$Limbic
length(which(rownames(cort.ahba.limbic) == ahba.probes.sort$gene_symbol))
rownames(cort.ahba.limbic) <- as.character(ahba.probes.sort$entrez_id)
cort.ahba.limbic.sort <- cort.ahba.limbic[order(rownames(cort.ahba.limbic)),]

cort.ahba.limbic.subsetByBernard <- cort.ahba.limbic.sort[rownames(cort.ahba.limbic.sort) %in% rownames(bernard.limbic.sort),]
ahba.probes.sort.subsetByBernard <- ahba.probes.sort[rownames(cort.ahba.limbic.sort) %in% rownames(bernard.limbic.sort),]
bernard.limbic.subsetByAHBA      <- bernard.limbic.sort[rownames(bernard.limbic.sort) %in% rownames(cort.ahba.limbic.sort),]
length(which(rownames(bernard.limbic.subsetByAHBA) == rownames(cort.ahba.limbic.subsetByBernard)))

# correlation of AHBA limbic cortex to BERNARD limbic cortex
cor.test(bernard.limbic.subsetByAHBA$logFC, cort.ahba.limbic.subsetByBernard$logFC)

# Plot 
# ------------------------
ahba.limbic.cort.entrez <- rownames(cort.ahba.limbic)[intersect(which(cort.ahba.limbic$logFC > 0), which(cort.ahba.limbic$adj.P.Val <=.01))]
monk.replicationPlot(ahba.dat = cort.ahba.limbic.subsetByBernard, 
                     ahba.sig.genes = as.character(ahba.limbic.cort.entrez),
                     rep.dat = bernard.limbic.subsetByAHBA,
                     xlabel = 'AHBA Limbic Cortex Log2 Fold Change (n=6)',
                     ylabel = 'Bernard Limbic Cortex Log2 Fold Change (n=6)',
                     xlimit = 2.5, 
                     ylimit = 2.5)


bspan.limbic.sort     <- limbic.stats[order(rownames(limbic.stats)),]
gene.rows.sort        <- ahba.probes.sort[order(as.character(ahba.probes.sort$gene_symbol)),]
gene.rows.sort.subset <- gene.rows.sort[gene.rows.sort$gene_symbol %in% rownames(bspan.limbic.sort),]
length(which(gene.rows.sort.subset$gene_symbol %in% rownames(bspan.limbic.sort)))
rownames(bspan.limbic.sort) <- gene.rows.sort.subset$entrez_id


bspan.limbic.subsetByBernard <- bspan.limbic.sort[rownames(bspan.limbic.sort) %in% rownames(bernard.limbic.sort),]
bspan.limbic.subsetByBernard <- bspan.limbic.subsetByBernard[order(as.character(rownames(bspan.limbic.subsetByBernard))),]
bernard.limbic.subsetByBSPAN <- bernard.limbic.sort[rownames(bernard.limbic.sort) %in% rownames(bspan.limbic.subsetByBernard),]
bernard.limbic.subsetByBSPAN <- bernard.limbic.subsetByBSPAN[order(as.character(rownames(bernard.limbic.subsetByBSPAN))),]
# correlation of BRAINSPAN limbic cortex to BERNARD limbic cortex
cor.test(bernard.limbic.subsetByBSPAN$logFC, bspan.limbic.subsetByBernard$logFC)

bspan.limbic.cort.entrez <- rownames(bspan.limbic.sort)[intersect(which(bspan.limbic.sort$logFC > 0), which(bspan.limbic.sort$adj.P.Val <=.01))]
# Plot 
# ----------------------
monk.replicationPlot(ahba.dat = bspan.limbic.subsetByBernard, 
                     ahba.sig.genes = as.character(bspan.limbic.cort.entrez),
                     rep.dat = bernard.limbic.subsetByBSPAN,
                     xlabel = 'BSPAN Limbic Striatum Log2 Fold Change (n=6)',
                     ylabel = 'Bernard Limbic Cortex Log2 Fold Change (n=111)',
                     xlimit = 2.5, 
                     ylimit = 2.5)


# Convert the primate probe names to AHBA gene names
# -------------------------
avg.bernard.expr.plot <- averaged.data.unique #[which(rownames(averaged.data.unique) %in% name_arr),]

# Plot gene expression 
# -------------------------
head(averaged.data.unique)
bernard.plot.regions <- c("Orbitofrontal cortex", "Anterior cingulate gyrus", "Temporal area", "Middle temporal area", "Dorsolateral Prefrontal Cortex",
                          "Primary auditory ctx", "Primary somatosensory cortex", "Primary motor cortex", 'V1', 'V2')
plot.bernard(averaged.data.unique, colnames(averaged.data.unique), '6751', bernard.plot.regions, ylimit=2)
# -------------------------


# Determine which limbic associated genes gene overlap between AHBA human/primates
# -------------------------
# Primate cortical 
bernard.significant.limbic.cortex.idxs <- intersect(which(bernard.limbic.subsetByAHBA$logFC > 0), which(bernard.limbic.subsetByAHBA$adj.P.Val <= .01))
bernard.sig.limbic.genes               <- rownames(bernard.limbic.subsetByAHBA)[bernard.significant.limbic.cortex.idxs]
# Primate striatal 
blueprint.sig.limbic.striatum <- intersect(which(blueprint.nacc.ahba.subset$logFC > 0), which(blueprint.nacc.ahba.subset$adj.P.Val <= .01))
blueprint.sig.limbic.genes    <- rownames(blueprint.nacc.ahba.subset)[blueprint.sig.limbic.striatum]
# primate limbic cortico-striatal 
monkey.sig.limbic.genes <- intersect(blueprint.sig.limbic.genes, bernard.sig.limbic.genes)
ahba.sig.limbic.genes   <- limbic.entrez

# Hypergeometric test, primate limbic network
# -------------------------
num.overlap        <- length(intersect(blueprint.sig.limbic.genes, bernard.sig.limbic.genes))
num.genes          <- dim(bernard.limbic.table.subset)[1]
white.genes        <- length(blueprint.sig.limbic.genes)
black.genes        <- num.genes-length(bernard.sig.limbic.genes)
genes.drawn        <- length(bernard.sig.limbic.genes)
print(num.overlap)
phyper(num.overlap-1, white.genes, black.genes, genes.drawn, lower.tail = F)

# Hypergeometric test, primate limbic network to AHBA limbic network
# -------------------------
length(intersect(ahba.sig.limbic.genes, monkey.sig.limbic.genes))
num.overlap        <- length(intersect(ahba.sig.limbic.genes, monkey.sig.limbic.genes))
num.genes          <- dim(bernard.limbic.subsetByAHBA)[1]
white.genes        <- length(monkey.sig.limbic.genes)
black.genes        <- num.genes-length(monkey.sig.limbic.genes)
genes.drawn        <- length(ahba.sig.limbic.genes)
phyper(num.overlap-1, white.genes, black.genes, genes.drawn, lower.tail = F)


# Hypergeometric test, primate limbic network to replication limbic network
# -------------------------
replication.limbic.entrez <- ahba.probes$entrez_id[ahba.probes$gene_symbol %in% replication.limbic.genes]
num.overlap        <- length(intersect(replication.limbic.entrez, monkey.sig.limbic.genes))
num.genes          <- dim(bernard.limbic.subsetByAHBA)[1]
white.genes        <- length(monkey.sig.limbic.genes)
black.genes        <- num.genes-length(monkey.sig.limbic.genes)
genes.drawn        <- length(replication.limbic.genes)
phyper(num.overlap-1, white.genes, black.genes, genes.drawn, lower.tail = F)

# Somatomotor Differential expression in macaque cortex
# ------------------------------
cur_contrast <- paste('-0.1428571*Orbitofrontal_cortex-0.1428571*Anterior_cingulate_gyrus-0.1428571*Middle_temporal_area-',
                      '0.1428571*Dorsolateral_Prefrontal_Cortex+0.33333*Primary_auditory_ctx',
                      '+0.33333*Primary_somatosensory_cortex-0.1428571*V1-0.1428571*V2-0.1428571*Temporal_area+0.33333*Primary_motor_cortex', sep = '')
colnames(design) <- gsub(' ', '_', gsub('fac', '', colnames(design)))
corfit           <- duplicateCorrelation(averaged.data.unique, design, block=donor.arr)
cmtx         <- makeContrasts(cur_contrast, levels=design ) # Make a contrast matrix
tmplm        <- lmFit( averaged.data.unique, design, block=donor.arr, correlation=corfit$consensus.correlation ) # Fit the model to the data
bernard.sommot.fit.cort <- eBayes( contrasts.fit( tmplm, cmtx ) )
bernard.sommot.table    <- topTable(bernard.sommot.fit.cort, num=Inf)
bernard.sommot.stats    <- bernard.sommot.table[order(rownames(bernard.sommot.table)),]

write.csv(file=paste0(base.dir, '/output_files/bernard_cort_sommot.csv'), x=bernard.sommot.stats)

bernard.sommot.stats[rownames(bernard.sommot.stats) =='5816',] # pvalb
bernard.limbic.table[rownames(bernard.limbic.table) =='6750',] # sst
bernard.limbic.table[rownames(bernard.limbic.table) =='6751',] # sstr1

                          
# compare blueprint to AHBA
rownames(bernard.sommot.stats) <- as.character(rownames(bernard.sommot.stats))
bernard.sommot.sort  <- bernard.sommot.stats[order(rownames(bernard.sommot.stats)),]

ahba.probes      <- all_data[[1]]$probes_collapse
ahba.probes.sort <- ahba.probes[order(ahba.probes$gene_symbol),]
cort.ahba.sommot <- cort.foldchanges.n6$stats$SomMot
length(which(rownames(cort.ahba.sommot) == ahba.probes.sort$gene_symbol))
rownames(cort.ahba.sommot) <- as.character(ahba.probes.sort$entrez_id)
cort.ahba.sommot.sort      <- cort.ahba.sommot[order(rownames(cort.ahba.sommot)),]

cort.ahba.sommot.subsetByBernard <- cort.ahba.sommot.sort[rownames(cort.ahba.sommot.sort) %in% rownames(bernard.sommot.sort),]
ahba.probes.sort.subsetByBernard <- ahba.probes.sort[rownames(cort.ahba.sommot.sort) %in% rownames(bernard.sommot.sort),]
bernard.sommot.subsetByAHBA      <- bernard.sommot.sort[rownames(bernard.sommot.sort) %in% rownames(cort.ahba.sommot.sort),]
length(which(rownames(bernard.sommot.subsetByAHBA) == rownames(cort.ahba.sommot.subsetByBernard)))

# correlation of AHBA sommot cortex to BERNARD sommot cortex
cor.test(bernard.sommot.subsetByAHBA$logFC, cort.ahba.sommot.subsetByBernard$logFC)


# compare blueprint to AHBA
bspan.sommot.sort     <- sommot.stats[order(rownames(sommot.stats)),]
gene.rows.sort        <- ahba.probes.sort[order(as.character(ahba.probes.sort$gene_symbol)),]
gene.rows.sort.subset <- gene.rows.sort[gene.rows.sort$gene_symbol %in% rownames(bspan.sommot.sort),]
length(which(gene.rows.sort.subset$gene_symbol %in% rownames(bspan.sommot.sort)))
rownames(bspan.sommot.sort) <- gene.rows.sort.subset$entrez_id


bspan.sommot.subsetByBernard <- bspan.sommot.sort[rownames(bspan.sommot.sort) %in% rownames(bernard.sommot.sort),]
bspan.sommot.subsetByBernard <- bspan.sommot.subsetByBernard[order(as.character(rownames(bspan.sommot.subsetByBernard))),]
bernard.sommot.subsetByBSPAN <- bernard.sommot.sort[rownames(bernard.sommot.sort) %in% rownames(bspan.sommot.subsetByBernard),]
bernard.sommot.subsetByBSPAN <- bernard.sommot.subsetByBSPAN[order(as.character(rownames(bernard.sommot.subsetByBSPAN))),]
# correlation of BRAINSPAN sommot cortex to BERNARD sommot cortex
cor.test(bernard.sommot.subsetByBSPAN$logFC, bspan.sommot.subsetByBernard$logFC)

# Plot 
# ------------------------
bernard.sommot.stats.subset           <- bernard.sommot.stats[which(rownames(bernard.sommot.stats) %in% name_arr),]
rownames(bernard.sommot.stats.subset) <- gene.arr
write.csv(x=bernard.sommot.stats.subset[order(rownames(bernard.sommot.stats.subset)),], file=paste(base.dir, '/output_files/Bernard_somMot_Cort_DEX.csv', sep=''))






# --------------------------------------------------
# --------------------------------------------------
# Determine layer-specific expression of limbic network genes
# --------------------------------------------------
# --------------------------------------------------

# Run layer specific analysis for genes that are expressed in limbic cortico-striatal loop in AHBA+replication data
# -------------------------
bernard.col.names <- colnames(bernard.expr.unique)
layer.array <- NULL
subj.arr    <- NULL
# extract layer and subject info from column names
for (col.name in bernard.col.names){
  layer       <- strsplit(strsplit(col.name, ',')[[1]][2], ';')[[1]][1]
  layer.name  <- gsub(' ','',layer)
  layer.array <- c(layer.array, layer.name)
  subj        <- strsplit(strsplit(col.name, ',')[[1]][2], ';')[[1]][2]
  subj.name   <- gsub(' ','',subj)
  subj.arr    <- c(subj.arr, subj.name)
}
# Only examine OFC and ACC regions
ofc.acc.idxs   <- grep('Orbitofrontal|Anterior cingulate', bernard.col.names)
ofc.acc.layers <- layer.array[ofc.acc.idxs]
ofc.acc.eset   <- bernard.expr.unique[,ofc.acc.idxs]
ofc.acc.subj   <- subj.arr[ofc.acc.idxs]


# Subset significant AHBA limbic network genes, based on Bakken mapping
# -------------------------
ahba.limbic.entrez <- all_data[[1]]$probes_collapse$entrez_id[which(all_data[[1]]$probes_collapse$gene_symbol %in% ahba.limbic.genes)]
entrez.arr <- NULL
ct <- 1
for (affy in rownames(ofc.acc.eset)){
  ct <- ct + 1
  print(ct)
  cur.entrez <- bakken.gene.mapping[which(bakken.gene.mapping$probeid == affy),]$human_entrezid
  entrez.arr <- c(entrez.arr, cur.entrez)
}
rownames(ofc.acc.eset) <- entrez.arr
# -------------------------
off.acc.limbic.eset <- ofc.acc.eset[which(rownames(ofc.acc.eset) %in% limbic.entrez),]
ofc.genes <- NULL
for (entrez in rownames(off.acc.limbic.eset)){
  ent <- as.numeric(entrez)
  cur.gene <- as.character(all_data[[1]]$probes_collapse$gene_symbol[which(all_data[[1]]$probes_collapse$entrez_id == ent)])
  ofc.genes <- c(ofc.genes, cur.gene)
}
rownames(off.acc.limbic.eset) <- ofc.genes
# -------------------------

# Function to mean-normalize the genes across layers. 
# -------------------------
meanNormalize <- function(x) { 
  out.x <- x - mean(x) 
  return(out.x)
}
# -------------------------
# group samples by layer
# -------------------------
fac              <- as.factor(ofc.acc.layers)
design           <- model.matrix(~fac + 0)
colnames(design) <- gsub(' ', '_', gsub('fac', '', colnames(design)))
corfit           <- duplicateCorrelation(off.acc.limbic.eset, design, block=ofc.acc.subj)
# -------------------------
cur.contrast <- paste('1*', choi, '-', mult.term, '*', paste(o.chois, collapse = paste('-', mult.term, '*', sep = '')), sep = '')
cur.contrast <- '.5*layer2+.5*layer3-.5*layer5-.5*layer6'
# Make the contrast matrix
# -------------------------
cmtx         <- makeContrasts(cur.contrast, levels = design) 
# -------------------------
tmplm            <- lmFit( off.acc.limbic.eset, design, block=ofc.acc.subj, correlation=corfit$consensus.correlation) # Fit the model to the data
bernard.fit.cort <- eBayes( contrasts.fit( tmplm, cmtx ) )
dex.table <- topTable(bernard.fit.cort, number = Inf, adjust = "BH", sort.by = "none")
write.csv(x=dex.table[order(rownames(dex.table)),], file=paste(base.dir, '/output_files/Layer_specific_limbic.csv', sep=''))
# -------------------------


# mean-normalize the data for plotting
# ---------------------
ofc.acc.limbic.eset.meanNorm <- apply(off.acc.limbic.eset, 1, meanNormalize)
plot.mat <- ofc.acc.limbic.eset.meanNorm[,which(dex.table$adj.P.Val < .05)] # genes which are differentially expressed between deep/superficial layers
# ---------------------
# Color-scale for the plot
# ---------------------
mycol_red  <- colorpanel(n=100, low="white", mid="orange", high="red")
mycol_blue <- colorpanel(n=100, low="blue", mid="turquoise", high="white")
mycol      <- c(mycol_blue, mycol_red)
color_breaks = c(seq(-2.5,2.5,length=length(mycol)+1))
# ---------------------
heatmap.2(as.matrix(plot.mat), col = mycol, trace='none', 
          distfun = function(x) dist(x,method = 'euclidean'),
          breaks=color_breaks, symm=F,symkey=F,symbreaks=T, scale="none")



# ---------------------------------
# ---------------------------------
# Cell-specific enrichments, using data from Doyle et al. 2008
# Supplementary Figure 5
# ---------------------------------
# ---------------------------------

# -------------------
# Read the series matrix file
# -------------------
dat.in      <- getGEO(filename=paste(base.dir, 'data/GSE13379/GSE13379_series_matrix.txt', sep = ''), GSEMatrix=TRUE)
save(file=paste(base.dir, 'data/GSE13379/GSE13379_series_matrix.Rdata', sep = ''), x=dat.in)
#load(file=paste(base.dir, 'data/GSE13379/GSE13379_series_matrix.txt', sep = ''))
cell_source     <- dat.in@phenoData@data$source_name_ch1
mouse.cell.data <- log2(dat.in@assayData$exprs)


# -------------------
# map mouse Affymetrix IDs to ensemble IDs
# -------------------
mouse              <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
attributes         <- c("affy_mouse430_2","ensembl_gene_id",'hgnc_id')
mouse.affy2ens     <- getBM(attributes = attributes, values = TRUE, mart = mouse, bmHeader=FALSE)
use.mouse.affy2ens <- mouse.affy2ens[mouse.affy2ens$affy_mouse430_2 != '',]
# -------------------

# -------------------
# Get Mouse homologs of human genes
# -------------------
human      <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- c("ensembl_gene_id","mmusculus_homolog_ensembl_gene","mmusculus_homolog_perc_id_r1")
attributes <- c(attributes,"mmusculus_homolog_orthology_type","mmusculus_homolog_subtype", "mmusculus_homolog_associated_gene_name")
orth.mouse <- getBM( attributes, filters="with_mmusculus_homolog", values = TRUE, mart = human, bmHeader=FALSE)
# -------------------
# Map ensemble IDs to valid entrez IDs
# -------------------
entrez2ensemble <- getBM(attributes = c('entrezgene', 'hgnc_symbol', 'entrezgene_trans_name','ensembl_gene_id'), values = TRUE, mart = human, bmHeader=FALSE)
entrez.dict     <- entrez2ensemble[!is.na(entrez2ensemble$entrezgene),]
use.entrezs     <- all_data$`9861`$probes_collapse$entrez_id[all_data$`9861`$probes_collapse$gene_symbol %in% ahba.limbic.genes]
# -------------------

# Genes that have valid ensembl ids
# dont assume that the HGNC symbol is uniform across AHBA and the 'entrez.dict'
# -------------------
matching.genes      <- entrez.dict[which(entrez.dict$entrezgene %in% limbic.entrez),] # genes from the ahba limbic list that are in the entrez dictionary
nonmatching.entrezs <- use.entrezs[!use.entrezs %in% entrez.dict$entrezgene] # genes that don't
# For each gene in the entrez dictionary, get the corresponding gene name in the AHBA data
ahba.gene.arr <- NULL
for (entrez in matching.genes$entrezgene){
  cur.gene      <- as.character(all_data$`9861`$probes_collapse$gene_symbol[all_data$`9861`$probes_collapse$entrez_id == entrez])
  ahba.gene.arr <- c(ahba.gene.arr, cur.gene)
}
# "matching.genes" contains HGNC/ENTREZ/ENSEMBL/AHBA Gene Name
matching.genes <- cbind(matching.genes, ahba.gene.arr)


# Now match the limbic associated genes to the mouse homolog list based on Ensembl IDs 
# ---------------
mouse.ensembles <- orth.mouse[orth.mouse$ensembl_gene_id %in% matching.genes$ensembl_gene_id,] # genes with a valid human/mouse ensembl correspondence 
mouse.genes     <- NULL
# get the corresponding mouse gene names for the mouse ensembl ids
for (ens in mouse.ensembles$ensembl_gene_id){
  mouse.genes <- c(mouse.genes, as.character(matching.genes[matching.genes$ensembl_gene_id == ens,]$ahba.gene.arr))
}
mouse.ensembles <- cbind(mouse.genes, mouse.ensembles)


# Map mouse ensembl ids to affymetrix names
# ----------------------
sig.affy.probes <- use.mouse.affy2ens[use.mouse.affy2ens$ensembl_gene_id %in% mouse.ensembles$mmusculus_homolog_ensembl_gene,]
final.genes <- NULL
for (ens in sig.affy.probes$ensembl_gene_id){
  print(as.character(mouse.ensembles[mouse.ensembles$mmusculus_homolog_ensembl_gene == ens,]$mouse.genes))
  final.genes <- c(final.genes, as.character(mouse.ensembles[mouse.ensembles$mmusculus_homolog_ensembl_gene == ens,]$mouse.genes[1]))
}
sig.affy.probes <- cbind(final.genes, sig.affy.probes)

# put the expression info in the correct order for input to the collapseRows function
# ----------------------
sig.mouse.eset <- NULL
for (affy in sig.affy.probes$affy_mouse430_2 ){
  cur.row <- mouse.cell.data[rownames(mouse.cell.data) == affy,]
  sig.mouse.eset <- rbind(sig.mouse.eset, cur.row)
}
rownames(sig.mouse.eset) <- 1:dim(sig.mouse.eset)[1]
# ----------------------
# For probes with non-unique matches to human genes, select the one with the max mean
# ----------------------
out <- collapseRows(sig.mouse.eset, rowGroup = as.character(sig.affy.probes$final.genes), rowID = rownames(sig.mouse.eset),
                    method = 'MaxMean', connectivityBasedCollapsing = FALSE)
sig.mouse.orthologs           <- out$datETcollapsed
colnames(sig.mouse.orthologs) <- as.character(cell_source)
# ----------------------

# Only use TRAP method data
# ----------------------
use.mouse.eset <- t(sig.mouse.orthologs[,grep('TRAP', colnames(sig.mouse.orthologs))])
process_rows <- function(string){
  split_string <- strsplit(string, ',')[[1]]
  new_string   <- paste(split_string[1:length(split_string)-1], collapse='')
  return(new_string)
}
new_rownames <- unlist(lapply(rownames(use.mouse.eset), process_rows))
# ----------------------
# Average expression data from the same cell type
# ----------------------
out <- collapseRows(use.mouse.eset, rowGroup = new_rownames, rowID = rownames(use.mouse.eset),
                    method = 'Average', connectivityBasedCollapsing = FALSE)
use.mouse.eset.collapse <- out$datETcollapsed
# ----------------------

# WGCNA analysis
# ----------------------
disableWGCNAThreads()
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft    <- pickSoftThreshold(as.matrix(use.mouse.eset.collapse), powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


# Select the softpower level
# -----------------
softPower     <- 7;
adjacency.mat <- adjacency(use.mouse.eset.collapse, power = softPower);


# Turn adjacency into topological overlap
# --------------------------------
TOM           <- TOMsimilarity(adjacency.mat);
colnames(TOM) <- rownames(adjacency.mat)
rownames(TOM) <- rownames(adjacency.mat)
dissTOM       <- 1-TOM
# Call the hierarchical clustering function
geneTree      <- hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# set minimum module size to 10
# --------------------------------
minModuleSize = 10;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


# Convert numeric lables into colors
# --------------------------------
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
# Calculate eigengenes
MEList        <- moduleEigengenes(use.mouse.eset.collapse, colors = dynamicColors)
MEs           <- MEList$eigengenes
rownames(MEs) <- rownames(use.mouse.eset.collapse)
# Calculate dissimilarity of module eigengenes
MEDiss        <- 1-cor(MEs);
# Cluster module eigengenes
METree        <- hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")


# Merge simillar modules 
# ------------------
MEDissThres <- 0.3
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge        <- mergeCloseModules(use.mouse.eset.collapse, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors <- merge$colors;
# Eigengenes of the new merged modules:
mergedMEs    <- merge$newMEs;
rownames(mergedMEs) <- rownames(use.mouse.eset.collapse)


# Rename to moduleColors
moduleColors <- mergedColors
# Construct numerical labels corresponding to the colors
colorOrder   <- c("grey", standardColors(50));
moduleLabels <- match(moduleColors, colorOrder)-1;
MEs <- mergedMEs;
# Save module colors and labels for use in subsequent parts


# Get the data ready for plotting 
# -----------------
cell.plot.order <- c("Cort+ Neurons Cortex", "Layer 5a Neurons Cortex", "Layer 5b Neurons Cortex", "Layer 6 Neurons Cortex",
                     "Cck+ Neurons Cortex", "Pnoc+ Neurons Cortex", "Drd1+ Medium Spiny Neurons Striatum", "Drd2+ Medium Spiny Neurons Striatum",
                     "Cholinergic Neurons Corpus Striatum", "Cholinergic Neurons Basal Forebrain", "Cholinergic Neurons Spinal Cord",
                     "Motor Neurons Brainstem", "Granule Neurons Cerebellum", "Unipolar Brush Neurons Cerebellum", "Golgi Neurons Cerebellum", "Stellate and Basket Neurons Cerebellum",
                     "Purkinje Neuron Cerebellum", "Mature Oligodendrocytes Cortex", "Mixed Oligodendroglia Cortex",
                     "Mature Oligodendrocytes Cerebellum", "Mixed Oligodendroglia Cerebellum",
                     "Astrocytes Cerebellum", "Astrocytes Cortex", "Bergmann Glia Cerebellum")
cell.colors <- rownames(mergedMEs)
cell.colors[intersect(grep('Cortex',cell.colors), grep('Neurons',cell.colors))] <- 'dodgerblue4'
cell.colors[intersect(grep('Striatum',cell.colors), grep('Neurons',cell.colors))] <- 'dodgerblue1'
cell.colors[intersect(grep('Brainstem',cell.colors), grep('Neurons',cell.colors))] <- 'deepskyblue'
cell.colors[intersect(grep('Cerebellum',cell.colors), grep('Neurons',cell.colors))] <- 'deepskyblue'
cell.colors[intersect(grep('Cerebellum',cell.colors), grep('Neuron',cell.colors))] <- 'deepskyblue'
cell.colors[intersect(grep('Basal Forebrain',cell.colors), grep('Neurons',cell.colors))] <- 'deepskyblue'
cell.colors[intersect(grep('Spinal Cord',cell.colors), grep('Neurons',cell.colors))] <- 'deepskyblue'
cell.colors[grep('Oligodendrocytes',cell.colors)] <- 'slateblue1'
cell.colors[grep('Oligodendroglia',cell.colors)] <- 'slateblue1'
cell.colors[grep('Astrocytes',cell.colors)] <- 'orchid1'
cell.colors[grep('Glia',cell.colors)] <- 'orchid1'


# Create a plot for each module 
# ----------------------
color.scale <- factor(c('dodgerblue4', 'dodgerblue1', 'deepskyblue', 'slateblue1', 'orchid1'), levels = c('dodgerblue4', 'dodgerblue1', 'deepskyblue', 'slateblue1', 'orchid1'))
color.scale <- c('dodgerblue4', 'dodgerblue1', 'deepskyblue', 'slateblue1', 'orchid1')
for (mod in names(mergedMEs)){
  plot.me <- as.data.frame(cbind(mergedMEs[[mod]], rownames(mergedMEs), cell.colors))
  plot.me$V1 <- as.numeric(as.character(plot.me$V1))
  plot.me$V2 <- factor(plot.me$V2, levels = cell.plot.order)
  plot.me$cell.colors <- factor(plot.me$cell.colors, levels = color.scale)
  
  max.axis = round(max(abs(plot.me$V1)),1) +.1
  ggplot(plot.me, aes(V2, V1, fill=cell.colors)) +
    geom_bar(stat = "identity", position = "identity") + 
    scale_fill_manual(values=color.scale) + 
    theme_bw() + 
    scale_y_continuous(limits = c(-max.axis, max.axis)) + 
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  out.name <- paste(base.dir, 'figures/', mod, '.pdf', sep ='')
  ggsave(filename=out.name)
}




