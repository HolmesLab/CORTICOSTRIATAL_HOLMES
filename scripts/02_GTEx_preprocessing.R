library(DESeq2)
library(GenomicFeatures)
library(data.table)
library(edgeR)
library(biomartr) #



# Read raw, un-normalized gene read counts from GTEx
# --------------------------
base.dir <- '/Users/kevinanderson/PHD/PROJECTS/2017_CORTICOSTRIATAL_NATCOMM/'  # enter the directory containing the script repository
gtex.read <- as.data.frame(fread(paste0(base.dir, 'data/GTEx/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct'), header = TRUE))

# Read GTEx sample information
# --------------------------
gtex.sample.attrib <- fread(paste0(base.dir, 'data/GTEx/GTEx_v7_Annotations_SampleAttributesDS.txt'), header = TRUE)
gtex.pheno         <- read.csv(paste0(base.dir, 'data/GTEx/GTEx_v7_Annotations_SubjectPhenotypesDS.txt'), sep='\t')

# Identify striatal samples
# --------------------------
striatal.regions   <- c("Brain - Caudate (basal ganglia)", "Brain - Nucleus accumbens (basal ganglia)", "Brain - Putamen (basal ganglia)")
striatal.reg.short <- c('Caudate', 'Nucleus_accumbens', 'Putamen')
gtex.striat.attrib <- gtex.sample.attrib[which(gtex.sample.attrib$SMTSD %in% striatal.regions),]


# Determine which genes are in all 3 regions
# ----------------------
covariate.df <- NULL
count.df     <- NULL
for (reg.idx in 1:3){ 
  # Read fully normalized expression values that were used in the FastQTL GTEx analyses
  # these will be used to subset the unnormalized count values
  # --------------------------
  reg <- striatal.reg.short[reg.idx]
  print(reg)
  tissue.file <- paste0(base.dir, 'data/GTEx/GTEx_Analysis_v7_eQTL_expression_matrices/Brain_', reg,'_basal_ganglia.v7.normalized_expression.bed')
  tissue.dat  <- read.table(tissue.file, comment.char='', header=TRUE)
  tissue.expr <- tissue.dat[,grep('GTEX', colnames(tissue.dat))]
  rownames(tissue.expr) <- tissue.dat$gene_id
  count.df[[reg]] <- tissue.expr
}
# only keep genes that are present in each of the 3 striatal regions
# -----------------------------
all.genes <- NULL
for (i in names(count.df)){
  all.genes <- c(all.genes, rownames(count.df[[i]]))
}
gene.table <- table(all.genes) # create table of gene counts
genes.keep <- names(gene.table)[which(gene.table == 3)]



# Extract count data for samples passing GTEx QC thresholds and genes in all three striatal regions
covariate.df <- NULL
count.df     <- matrix()
for (reg.idx in 1:3){ 
  # Read fully normalized expression values that were used in the FastQTL GTEx analyses
  # these will be used to subset the unnormalized count values
  # --------------------------
  reg <- striatal.reg.short[reg.idx]
  print(reg)
  tissue.file <- paste0(base.dir, 'data/GTEx/GTEx_Analysis_v7_eQTL_expression_matrices/Brain_', reg,'_basal_ganglia.v7.normalized_expression.bed')
  tissue.dat  <- read.table(tissue.file, comment.char='', header=TRUE)
  tissue.expr <- tissue.dat[,grep('GTEX', colnames(tissue.dat))]
  rownames(tissue.expr) <- tissue.dat$gene_id
  
  # Subset Count data to only include samples from the current region
  # --------------------------
  region          <- striatal.regions[reg.idx]
  region.sampids  <- gtex.sample.attrib$SAMPID[which(gtex.sample.attrib$SMTSD %in% region)]
  region.samp.att <- gtex.sample.attrib[which(gtex.sample.attrib$SMTSD %in% region),]
  
  # pull the reads for this region (from the 'gtex.read' variable)
  # --------------------------
  region.counts  <- gtex.read[,which(colnames(gtex.read) %in% region.sampids)]
  rownames(region.counts) <- gtex.read$Name
  
  # identify the genes that passed GTEx's QC pipeline
  # --------------------------
  tpm.use.genes    <- gtex.read$Name[which(gtex.read$Name %in% tissue.dat$gene_id)]
  # identify samples that passed GTEx's QC pipeline
  # --------------------------
  gtex.use.samples <- colnames(tissue.dat)[grep('GTEX', colnames(tissue.dat))] # gtex processed
  gtex.read.refs   <- gsub('-', '.', colnames(region.counts)) # count data
  use.samp.idxs    <- unlist(lapply(gtex.use.samples, grep, x=gtex.read.refs)) # idxs in the count data for usable samples
  
  
  # Subset count data so that they reflects GTEx QC methods
  # --------------------------
  # genes that pass GTEx QC and are present in all 3 striatal samples
  use.genes      <- intersect(tpm.use.genes, genes.keep)
  USE.reg.counts <- region.counts[which(rownames(region.counts) %in% use.genes), use.samp.idxs]
  
  # append the samples as columns in the data frame
  count.df <- cbind(count.df, as.data.frame(USE.reg.counts))

  # Get region-wise covariates
  # --------------------------
  covar.file        <- paste0(base.dir, 'data/GTEx/GTEx_Analysis_v7_eQTL_covariates/Brain_', reg,'_basal_ganglia.v7.covariates.txt')
  tissue.covars.tmp <- read.table(covar.file, comment.char='', header=TRUE)  
  tissue.covar      <- tissue.covars.tmp[,grep('GTEX', colnames(tissue.covars.tmp))]
  rownames(tissue.covar) <- tissue.covars.tmp$ID
  region <- rep(reg, dim(tissue.covar)[2])
  
  # phenotype information for the current tissue
  # --------------------------
  region.pheno <- gtex.pheno[gsub('-','.',gtex.pheno$SUBJID) %in% colnames(tissue.covar),]
  
  
  # Check the order of each data frame: count / covariate / phenotype
  # (make sure subject order is maintained before concatenating into a larger dataframe)
  # --------------------------
  print(length(which(gsub('-','.',region.pheno$SUBJID) == colnames(tissue.covar))))
  count.sub.names <- unlist(lapply(colnames(USE.reg.counts), function(x) paste(strsplit(x, '-')[[1]][1:2], collapse='.') ))
  length(which(count.sub.names == colnames(tissue.covar)))
  
  # transpose covariates
  t.tissue.covar <- as.data.frame(t(tissue.covar))
  
  # append regional covariates to the overall covirate dataframe
  #sex      <- t.tissue.covar$sex
  #platform <- t.tissue.covar$platform
  sub_num      <- rownames(t.tissue.covar)
  covar.out    <- cbind(t.tissue.covar, sub_num, region, region.pheno)
  covariate.df <- rbind(covariate.df, covar.out)
}
count.df$count.df <- NULL # get rid of column header with NAs

# Create a gene name dictionary with entrez id/hgnc/ensembl
split_me <- function(x){ 
  tmp <- strsplit(as.character(x), '[.]')
  return(unlist(tmp)[1])
}
ensembl.base <- unlist(lapply(rownames(count.df), split_me))
rownames(count.df) <- ensembl.base

# -----------------------
# Dictionary with ensembl/entrez gene information
# -----------------------
library(biomaRt) #
mart        <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_names  <- getBM(filters="ensembl_gene_id",  attributes=c("hgnc_symbol", "ensembl_gene_id", "entrezgene", "description"),
                     values=ensembl.base, mart=mart)
# ------------------


# Convert ensembl gene ids to gene symbols/entrez ids
# ------------------
ensembl2gene        <- matrix(NA, length(ensembl.base), 1)
ensembl2entrez      <- matrix(NA, length(ensembl.base), 1)
ct <- 1
for (ens in ensembl.base){
  print(paste(ct,'/', length(ensembl.base), sep = '')) # in-line feedback
  ct <- ct + 1
  ens_idxs   <- which(ensembl.base == ens)
  entrez     <- gene_names$entrezgene[which(gene_names$ensembl_gene_id == ens)][1]
  ensembl2entrez[ens_idxs] <- entrez
  gene_acro  <- gene_names$hgnc_symbol[which(gene_names$ensembl_gene_id == ens)][1]
  ensembl2gene[ens_idxs]   <- gene_acro
}
rowIDs <- 1:length(ensembl.base)


# Remove genes that don't map to an entrez id
# ------------------
valid_idxs     <- which(!is.na(ensembl2entrez))
groupID        <- ensembl2entrez[valid_idxs]
rowID          <- rownames(count.df)[valid_idxs]
gtex.data      <- count.df[valid_idxs,]
out            <- collapseRows(gtex.data, groupID, rowID, method='MaxMean', connectivityBasedCollapsing = FALSE)
gtex.collapsed.count <- out$datETcollapsed

arr <- NULL
age.arr <- NULL
for (sub in as.character(unique(covariate.df$sub_num)) ){
  arr <- c(arr, covariate.df[which(covariate.df$sub_num == sub),]$sex[1])
  age.arr <- c(age.arr, as.character(covariate.df[which(covariate.df$sub_num == sub),]$AGE[1]))
}

library(edgeR)
dge.counts <- DGEList(counts=gtex.collapsed.count, 
                      samples=covariate.df,
                      lib.size=colSums(gtex.collapsed.count),
                      genes=rownames(gtex.collapsed.count))

cpm  <- cpm(dge.counts)
lcpm <- cpm(dge.counts, log=TRUE)
# table(rowSums(dge.counts$counts==0)==9)

keep.exprs <- which(rowSums(cpm>1)>=3)  # genes with at least 3 counts across all the data
dge.counts <- dge.counts[keep.exprs,, keep.lib.sizes=FALSE]



library(RColorBrewer)
nsamples <- ncol(dge.counts)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(dge.counts, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

# transform age data to meaningful ordinal variable
age.ordinal <- as.character(dge.counts$samples$AGE)
age.ordinal[dge.counts$samples$AGE == '20-29'] = 1
age.ordinal[dge.counts$samples$AGE == '30-39'] = 2
age.ordinal[dge.counts$samples$AGE == '40-49'] = 3
age.ordinal[dge.counts$samples$AGE == '50-59'] = 4
age.ordinal[dge.counts$samples$AGE == '60-69'] = 5
age.ordinal[dge.counts$samples$AGE == '70-79'] = 6
dge.counts$samples$Age_ordinal <- as.numeric(age.ordinal)

# normalisation using trimmed mean of M-values (TMM)
dge.counts <- calcNormFactors(dge.counts, method="TMM")

formula <- formula(paste0('~0+region + sex + platform + Age_ordinal + ', paste(paste0('InferredCov', 1:15), collapse=' + ')))
design  <- model.matrix(formula, data=dge.counts$samples)
colnames(design) <- gsub('region|sub_num','',colnames(design))

contr.matrix <- makeContrasts(
  NAcc_vs_CaudoPut = +1*Nucleus_accumbens-0.5*Caudate-0.5*Putamen, 
  levels = colnames(design))

# Voom requires information from duplicateCorrelation and thus must be run at least twice.
# voom
v <- voom(dge.counts, design, plot=TRUE)
vfit <- lmFit(v, design)
cfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(cfit)

# Second voom
corfit <- duplicateCorrelation(v, design, block=dge.counts$samples$sub_num)
v2     <- voom(dge.counts, design, plot=FALSE, block=dge.counts$samples$sub_num, correlation=corfit$consensus)
vfit.2 <- lmFit(v2, design, block=dge.counts$samples$sub_num, correlation=corfit$consensus)
cfit.2 <- contrasts.fit(vfit.2, contrasts=contr.matrix)
cfit.2 <- eBayes(cfit.2)

# Third voom
corfit <- duplicateCorrelation(v2, design, block=dge.counts$samples$sub_num)
v2     <- voom(dge.counts, design, plot=FALSE, block=dge.counts$samples$sub_num, correlation=corfit$consensus)
vfit.2 <- lmFit(v2, design, block=dge.counts$samples$sub_num, correlation=corfit$consensus)
cfit.2 <- contrasts.fit(vfit.2, contrasts=contr.matrix)
cfit.2 <- eBayes(cfit.2)

nacc.contrast <- topTable(cfit.2, number=Inf, sort.by="none")

# save output for later reading
save(file=paste(base.dir, 'data/GTEx/nacc.contrast.Rdata', sep = ''), x=nacc.contrast)
save(file=paste(base.dir, 'data/GTEx/dge.counts.Rdata', sep = ''), x=dge.counts)

