# Function library for corticostriatal analyses

# Written by: Kevin M. Anderson
# Contact:    kevin.anderson@yale.edu 

striatExprOfCortGeneSets <- function(cort.foldchanges.n6, striat.out.ind){
  all.striat.out <- NULL
  for (donor in donor.nums){
    cur.dat  <- as.data.frame(striat.out.ind[[donor]])
    expr.arr <- NULL
    for (reg in use.striat.networks){
      cur.genes   <- cort.foldchanges.n6$q01_genes[[reg]]
      within.mean <- mean(cur.dat[[reg]][rownames(cur.dat) %in% cur.genes])
      names(within.mean) <- reg
      o.regs      <- use.striat.networks[which(!use.striat.networks %in% reg)]
      
      betw.means <- NULL
      for (oreg in o.regs){
        betw.mean  <- mean(cur.dat[[oreg]][rownames(cur.dat) %in% cur.genes])
        betw.means <- c(betw.means, betw.mean)
      }
      names(betw.means) <- o.regs
      expr.arr <- c(within.mean, betw.means)
      group    <- c('within', rep('between', 1, 4))
      donor    <- rep(donor, 1, 5)
      network  <- rep(reg, 1, 5)
      dat.out  <- cbind(expr.arr, group, donor, network)
      colnames(dat.out) <- c('expr', 'group', 'donor', 'network')
      rownames(dat.out) <- NULL
      dat.out  <- as.data.frame(dat.out)
      dat.out$expr <- as.numeric(as.character(dat.out$expr))
      all.striat.out <- rbind(all.striat.out, dat.out)
    }
    # Expression of Som/Mot genes in VentAttn striatum
    reg <- 'VentAttn'
    cur.genes   <- cort.foldchanges.n6$q01_genes$SomMot
    within.mean <- mean(cur.dat[[reg]][rownames(cur.dat) %in% cur.genes])
    names(within.mean) <- reg
    o.regs      <- use.striat.networks[which(!use.striat.networks %in% reg)]
    betw.means  <- NULL
    for (oreg in o.regs){
      betw.mean  <- mean(cur.dat[[oreg]][rownames(cur.dat) %in% cur.genes])
      betw.means <- c(betw.means, betw.mean)
    }
    names(betw.means) <- o.regs
    expr.arr <- c(within.mean, betw.means)
    group    <- c('within', rep('between', 1, 4))
    donor    <- rep(donor, 1, 5)
    network  <- rep('VentAttn_SomMot', 1, 5)
    dat.out  <- cbind(expr.arr, group, donor, network)
    colnames(dat.out) <- c('expr', 'group', 'donor', 'network')
    rownames(dat.out) <- NULL
    dat.out  <- as.data.frame(dat.out)
    dat.out$expr <- as.numeric(as.character(dat.out$expr))
    all.striat.out <- rbind(all.striat.out, dat.out)
  }
  return(all.striat.out)
}



plotSuppFigure3 <- function(cort.foldchanges.n6, striat.foldchanges.n6){
  
  # Limbic Cortex - Limbic Striatum
  frame()
  setEPS()
  postscript(paste0(base.dir, 'figures/SuppFig3_LimbicCort_2_LimbicStriat.eps'))
  ols <- lm(cort.foldchanges.n6$fit_df$Limbic$coefficients ~ striat.foldchanges.n6$fit_df$Limbic$coefficients)
  xlimit=5
  ylimit=2
  print(ggplot() +
    geom_point( aes(x = striat.foldchanges.n6$fit_df$Limbic$coefficients, y = cort.foldchanges.n6$fit_df$Limbic$coefficients, color='black'), size = 2.5) +
    xlab('limbic striatum') + ylab('limbic cortex')  +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
    scale_x_continuous(limits = c(-1*xlimit, xlimit), breaks = number_ticks(11))  +
    scale_y_continuous(limits = c(-1*ylimit, ylimit), breaks = number_ticks(11))  +
    geom_abline(intercept = ols$coefficients[1],
                slope = ols$coefficients[2]))
  dev.off()
  # -------------------
  # -------------------
  # SomMot Cortex - Limbic Striatum
  frame()
  setEPS()
  postscript(paste0(base.dir, 'figures/SuppFig3_SomMotCort_2_LimbicStriat.eps'))
  ols <- lm(cort.foldchanges.n6$fit_df$SomMot$coefficients ~ striat.foldchanges.n6$fit_df$Limbic$coefficients)
  xlimit=5
  ylimit=2
  print(ggplot() +
    geom_point( aes(x = striat.foldchanges.n6$fit_df$Limbic$coefficients, y = cort.foldchanges.n6$fit_df$SomMot$coefficients, color='black'), size = 2.5) +
    xlab('limbic striatum') + ylab('SomMot cortex')  +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
    scale_x_continuous(limits = c(-1*xlimit, xlimit), breaks = number_ticks(11))  +
    scale_y_continuous(limits = c(-1*ylimit, ylimit), breaks = number_ticks(11))  +
    geom_abline(intercept = ols$coefficients[1],
                slope = ols$coefficients[2]))
  dev.off()
  # -------------------
  # -------------------
  # SomMot Cortex - SomMot Striatum
  setEPS()
  postscript(paste0(base.dir, 'figures/SuppFig3_SomMotCort_2_SomMotStriat.eps'))
  frame()
  ols <- lm(cort.foldchanges.n6$fit_df$SomMot$coefficients ~ striat.foldchanges.n6$fit_df$SomMot$coefficients)
  xlimit=4
  ylimit=2
  print(ggplot() +
    geom_point( aes(x = striat.foldchanges.n6$fit_df$SomMot$coefficients, y = cort.foldchanges.n6$fit_df$SomMot$coefficients, color='black'), size = 2.5) +
    xlab('SomMot striatum') + ylab('SomMot cortex')  +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
    scale_x_continuous(limits = c(-1*xlimit, xlimit), breaks = number_ticks(11))  +
    scale_y_continuous(limits = c(-1*ylimit, ylimit), breaks = number_ticks(11))  +
    geom_abline(intercept = ols$coefficients[1],
                slope = ols$coefficients[2]))
  dev.off()
  # -------------------
  # -------------------
  # Limbic Cortex - SomMot Striatum
  frame()
  setEPS()
  postscript(paste0(base.dir, 'figures/SuppFig3_LimbicCort_2_SomMotStriat.eps'))
  ols <- lm(cort.foldchanges.n6$fit_df$Limbic$coefficients ~ striat.foldchanges.n6$fit_df$SomMot$coefficients)
  xlimit=4
  ylimit=2
  print(ggplot() +
    geom_point( aes(x = striat.foldchanges.n6$fit_df$SomMot$coefficients, y = cort.foldchanges.n6$fit_df$Limbic$coefficients, color='black'), size = 2.5) +
    xlab('SomMot striatum') + ylab('Limbic cortex')  +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
    scale_x_continuous(limits = c(-1*xlimit, xlimit), breaks = number_ticks(11))  +
    scale_y_continuous(limits = c(-1*ylimit, ylimit), breaks = number_ticks(11))  +
    geom_abline(intercept = ols$coefficients[1],
                slope = ols$coefficients[2]))
  dev.off()
  # -------------------
}



ahba_diff_expression <- function(all_data, use.donors, use.regions, reg2yeo, dat2avg, rest.networks){
  
  out <- averageWithinCortNetworks(all_data=all_data, 
                                    use.regions=use.regions, 
                                    reg2yeo=reg2yeo, 
                                    use.donors=use.donors, # LH donors
                                    dat2avg=dat2avg, 
                                    regs=rest.networks)
  expr    <- out[[1]]
  regions <- out[[2]]
  donors  <- out[[3]]
  colnames(expr) <- paste(donors, regions, sep = '_')
  
  # Use limma to calculate differential expression for each network, relative to all others
  # ------------------------
  fac              <- as.factor(regions) # the network type of each column in 'expr'
  design           <- model.matrix(~0 + fac)
  colnames(design) <- gsub('fac', '', colnames(design))
  corfit           <- duplicateCorrelation(expr, design, block=donors) 
  # Use duplicate correlation to account for random effect of subject. The alternative is to use donor as an explicit fixed effect
  # in the model for a slight boost in statistical power. Random effects approach is used here for consistency with later GTEx/Brainspan 
  # analyses, where the random effects approach allows for more data retention (i.e. keep donors without complete pairwise data)
  
  # Calculate Differential expression for each region
  # ------------------------
  sig.genes       <- NULL
  cort.foldchange <- NULL
  cort.foldchange[['q01_genes']] <- list()
  for ( net in rest.networks ) {
    print(paste('Getting preferential expression for: ', net, sep = ''))
    
    # negative weight for the contrast matrix, depends on number of comparison networks
    mult.term    <- round(1/(length(region.names)-1),6)
    o.nets      <- region.names[region.names != net] # name of the other networks
    cur.contrast <- paste('1*', net, '-', mult.term, '*', paste(o.nets, collapse = paste('-', mult.term, '*', sep = '')), sep = '')
    
    # Make the contrast matrix
    cmtx <- NULL
    cmtx <- makeContrasts(contrasts=cur.contrast, levels=colnames(design))
    tmplm   <- lmFit(expr, design, block=donors, correlation=corfit$consensus.correlation ) # Fit the model to the data
    fit     <- eBayes(contrasts.fit( tmplm, cmtx ) )
    cort.foldchange[['fit_df']][[net]] <- fit
    tmp     <- topTable(fit, number=Inf)
    cort.foldchange[['stats']][[net]] <- tmp[order(rownames(tmp)),] 
    
    # Positive fold change, FDR corrected p<0.01
    pos.idxs       <- which(cort.foldchange$stats[[net]]$logFC > 0)
    adjusted.ps    <- which(cort.foldchange$stats[[net]]$adj.P.Val <= .01)
    genes.tmp      <- rownames(cort.foldchange$stats[[net]])[intersect(adjusted.ps, pos.idxs)]
    cort.foldchange[['q01_genes']][[net]] <- genes.tmp
    print(length(genes.tmp))
    sig.genes <- c(sig.genes, genes.tmp)
  }
  return(list(cort.foldchange, sig.genes))
}



avg_parcel_expression <- function(all_data, striatal.atlas, striatal.num, cortical.atlas, cortical.num, cort.use.regions){
  
  for ( donor in names(all_data)){
    print(paste('Averaging cort and striat data for: ', donor, sep = ''))
    
    # Get the striatal parcels that contain samples for this subject
    # ----------------------
    striat.regions.tmp <- unique(all_data[[donor]][[striatal.atlas]])
    striat.regions     <- sort(striat.regions.tmp[striat.regions.tmp != 0]) # Don't count areas that weren't assigned (i.e not 0)
    
    # Mean Normalized striatal expression values
    all_data[[donor]]$striat_expr <- averageStriatExpr(striat.region=striat.regions, 
                                                       data.struct=all_data[[donor]], 
                                                       choi.names=choi.names[[striatal.num]], 
                                                       samp.labels=all_data[[donor]][[striatal.atlas]], 
                                                       data.type='striat_meanNorm')
    # non-mean Normalized striatal expression values
    all_data[[donor]]$striat_expr_nonorm <- averageStriatExpr(striat.region=striat.regions, 
                                                              data.struct=all_data[[donor]], 
                                                              choi.names=choi.names[[striatal.num]], 
                                                              samp.labels=all_data[[donor]][[striatal.atlas]], 
                                                              data.type='all_striat_micros')
    
    # Mean Normalized cortical expression values
    all_data[[donor]]$cort_expr         <- averageCortExpr(use_regions=cort.use.regions, 
                                                           data_struct=all_data[[donor]], 
                                                           max_reg=cortical.num, 
                                                           sample_labels=all_data[[donor]][[cortical.atlas]], 
                                                           data_type='cort_meanNorm')
    # non-Mean Normalized cortical expression values
    all_data[[donor]]$cort_expr_nonorm  <- averageCortExpr(use_regions=cort.use.regions, 
                                                           data_struct=all_data[[donor]], 
                                                           max_reg=cortical.num, 
                                                           sample_labels=all_data[[donor]][[cortical.atlas]], 
                                                           data_type='all_cort_micros')
  }
  return(all_data)
}



# Cross-reference MNI coordinates with group atlas
# ----------------------------
cort_query <- function(data_struct, atlas, MNI_coords, rad, afni.dir) {
  
  # Initialize output arrays
  net_arr     <- NULL
  in_network  <- NULL
  out_network <- NULL
  # Iterate of all cortical sample locations
  # ---------------------
  for ( iter in 1:length(data_struct$cort_samples[,1]) ){
    use_coords       <- MNI_coords[MNI_coords$well_id %in% data_struct$cort_sample$well_id[iter],]
    coords           <- c(use_coords$mni_x, use_coords$mni_y, use_coords$mni_z)
    network          <- query_atlas(atlas, coords, 0, afni.dir)
    print(network)
    
    # Check neighboring voxels if needed
    if ( rad == 1 ){
      neigh_arr <- query_atlas(atlas, coords, 1, afni.dir)
      if ( length(neigh_arr[neigh_arr %in% 0]) == 27 ){
        out_network <- rbind(out_network, data_struct$cort_samples[iter,])
      } else if ( length(unique(neigh_arr[neigh_arr > 0])) == 1 ){
        if ( network == 0 ){
          use     <- neigh_arr[neigh_arr>0]
          network <- as.integer(names(sort(summary(as.factor(use))))[1])
        }
        in_network  <- rbind(in_network, data_struct$cort_samples[iter,]) 
      } else {
        network <- 999
        out_network <- rbind(out_network, data_struct$cort_samples[iter,])
      }
      net_arr   <- c(net_arr, network)
    } else {
      net_arr   <- c(net_arr, network)
    }
  }
  return(net_arr)
}



striat_query <- function(data_struct, atlas, atlas_conf, MNI_coords, afni.dir) {
  net_arr      <- NULL
  net_arr_conf <- NULL
  in_network   <- NULL
  conf_network <- NULL
  out_network  <- NULL
  
  for ( iter in 1:dim(data_struct$striat_samples)[1] ){
    print(iter)
    #coords       <- c(data_struct$striat_samples[iter,]$mni_x, data_struct$striat_samples[iter,]$mni_y, data_struct$striat_samples[iter,]$mni_z)
    use_coords <- MNI_coords[MNI_coords$well_id %in% data_struct$striat_samples$well_id[iter],]
    coords     <- c(use_coords$mni_x, use_coords$mni_y, use_coords$mni_z)
    
    network      <- query_atlas(atlas, coords, 0, afni.dir)
    network_conf <- query_atlas(atlas_conf, coords, 0, afni.dir)
    
    # Check neighboring voxels
    if ( network == 0 ){
      neigh_arr <- query_atlas(atlas, coords, 1, afni.dir)
      num_zeros <- length(neigh_arr[neigh_arr %in% 0])
      if ( num_zeros == 27 ){
        neigh_arr <- query_atlas(atlas, coords, 2, afni.dir)
        num_zeros <- length(neigh_arr[neigh_arr %in% 0])
        if ( num_zeros == 125 ){
          network <- 0 
        } else {
          use       <- neigh_arr[neigh_arr>0]
          use_table <- summary(as.factor(use))
          network   <- as.integer(names(sort(use_table)))
          network_conf <- 0
        }
      } else {
        use       <- neigh_arr[neigh_arr>0]
        use_table <- summary(as.factor(use))
        network   <- as.integer(names(sort(use_table)))
        network_conf <- 0
      }
    }
    print(paste(iter,':', network))
    net_arr      <- c(net_arr, network[1])
    net_arr_conf <- c(net_arr_conf, network_conf)
    
  }
  return(list(net_arr, net_arr_conf))
}




readAtlasOverlap <- function(all_data, filenames, atlas.dir){
  # Read in the atlas overlap information that aligns the MNI coordinates of each sample to the 
  # Choi/Yeo striatal and cortical atlases. These were calculated ahead of time in 'fyi_preprocess_data.R'
  # 
  # Args:
  #   all_data:  data structure containing microarray expr, as well as sample/ontology/probe info
  #   filenames: donor samples, corresponds to field names in 'all_data'
  #   atlas.dir: directory with all the atlas information
  #
  # Returns:
  #   Atlas overlap information for each subject
  
  
  for ( donor in filenames ){
    types  <- c('splitLabel_cort_', 'ChoiMNI152_striat_')
    n.regs <- c('7','17')
    
    for (type in types){
      for (nreg in n.regs){
        atlas.name <- paste0('/', type, donor, '_', nreg, 'net.csv')
        atlas.in   <- read.csv(file = paste(atlas.dir, atlas.name, sep=''))
        atlas.out  <- atlas.in$x
        cur.n      <- strsplit(atlas.name, '_')[[1]][1]
        use_name   <- paste(gsub('/', '', cur.n), nreg, sep = '_')
        print(use_name)
        all_data[[donor]][[use_name]] <- atlas.out
      }
    }
  }
  return(all_data)
}






calcParcelDist <- function(all_data, donors, regions){
  
  both.dists <- NULL
  for (donor in donors){
    donor.dists     <- matrix(NA, length(ordered.use.regions), length(ordered.use.regions))
    cort.samples    <- all_data[[donor]]$cort_samples
    cort.reg.labels <- all_data[[donor]]$splitLabel_17
    for (reg.i in 1:(length(ordered.use.regions)-1) ) {
      reg1        <- ordered.use.regions[reg.i]
      reg.idxs    <- which(cort.reg.labels %in% reg1)
      reg.samples <- cort.samples[reg.idxs,]
      avg.y.coord <- mean(reg.samples$mni_y)
      avg.x.coord <- mean(reg.samples$mni_x)
      avg.z.coord <- mean(reg.samples$mni_z)
      
      for ( reg.j in (reg.i + 1):length(ordered.use.regions) ){
        reg2         <- ordered.use.regions[reg.j]
        reg2.idxs    <- which(cort.reg.labels %in% reg2)
        reg2.samples <- cort.samples[reg2.idxs,]
        avg2.y.coord <- mean(reg2.samples$mni_y)
        avg2.x.coord <- mean(reg2.samples$mni_x)
        avg2.z.coord <- mean(reg2.samples$mni_z)
        
        
        # Euclidean distance between the two regions
        x <- (avg.x.coord - avg2.x.coord)^2
        y <- (avg.y.coord - avg2.y.coord)^2
        z <- (avg.z.coord - avg2.z.coord)^2
        dist <- sqrt(x + y + z)
        
        donor.dists[reg.i, reg.j] <- dist
        donor.dists[reg.j, reg.i] <- dist
      }
    }
    both.dists[[donor]] <- donor.dists
  }
  return(both.dists)
}





get_region_info <- function(all_data, filenames, reg_IDs, name){
  # Given a list of AHBA region ID numbers, subset the data for each subject
  
  # Args:
  #   all_data:  data structure containing microarray expression, as well as sample/ontology/probe info, 
  #   filenames: list of subject numbers, which correspond to field names of 'all_data'
  #   reg_IDs:   AHBA ontology IDs for the region being analyzed
  #   name:      string for distinguishing the output data
  #
  # Returns:
  #   Region specific versions of microarray/ ontology/ probe/ and sample information
  
  
  for ( donor in filenames ){
    print(paste('Getting ', name, ' data for: ', donor, sep = ''))
    
    # Get expression and sample information about cortical samples
    # -----------------
    idxs      <- all_data[[donor]]$raw_samp$structure_id %in% reg_IDs
    
    samp_name <- paste(name, '_samples', sep='')
    all_data[[donor]][[samp_name]]  <- all_data[[donor]]$raw_samp[idxs, ] # sample information
    
    acro_name <- paste(name, '_acros', sep='')
    all_data[[donor]][[acro_name]]  <- factor(all_data[[donor]][[samp_name]]$structure_acronym) # list of corresponding struct acronyms
    
    micro_name <- paste('all_', name, '_micros', sep='')
    all_data[[donor]][[micro_name]] <- as.data.frame(all_data[[donor]]$micro_collapsed)[, idxs]
    
    mean_name <- paste(name, '_meanNorm', sep='')
    all_data[[donor]][[mean_name]]  <- t(apply(all_data[[donor]][[micro_name]], 1, function(x) x-mean(x) ))
    
    pa_name <- paste('all_', name, '_pas', sep='')
    all_data[[donor]][[pa_name]] <- as.data.frame(all_data[[donor]]$pas_collapse)[, idxs]
    
  }
  return (all_data)
}


plot.bernard <- function(bernard.collapsed.data, blueprint.mature.cols, gene, regions2plot, ylimit){
  gene.idx <- which(rownames(bernard.collapsed.data) %in% gene)
  gene.dat <- bernard.collapsed.data[gene.idx,] 
  gene.dat <- gene.dat - mean(gene.dat)
  dat.out  <- NULL
  for (reg in regions2plot){
    reg.dat  <- gene.dat[grep(reg, names(gene.dat), ignore.case = FALSE)]
    reg.mean <- mean(reg.dat)
    reg.sd   <- sd(reg.dat)
    reg.se   <- sd(reg.dat) / sqrt(length(reg.dat))
    dat.out  <- rbind(dat.out, c(reg.mean, reg.sd, reg.se, reg))
  }
  dat.out <- as.data.frame(dat.out)
  rownames(dat.out) <- regions2plot
  dat.out$V1 <- as.numeric(as.character(dat.out$V1))
  dat.out$V2 <- as.numeric(as.character(dat.out$V2))
  dat.out$V3 <- as.numeric(as.character(dat.out$V3))
  colnames(dat.out) <- c('mean','sd','se','reg')
  dat.out$reg <- factor(dat.out$reg, levels = regions2plot)
  dat.out$mean <- dat.out$mean - mean(dat.out$mean)
  dodge <- position_dodge(width=0.9)
  limits <- aes(ymax = mean + se, ymin=mean - se)
  ggplot(dat.out, aes(x=reg, y=mean)) + geom_bar(position=dodge, stat="identity") + 
    geom_errorbar(limits, position=dodge, width=0.25) +
    ylim(-ylimit, ylimit) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    ggtitle(gene)
}



# Function to set the number of ticks on a ggplot graph
# --------------
number_ticks <- function(n) {function(limits) pretty(limits, n)}


# --------------
getNetworkSpecificExpr <- function( mat.dat, netGenes, region.names, gene.list ) {
  # Output structures
  dat.out   <- NULL
  aov.table <- NULL
  for ( choi in region.names ){
    # get the indices for the differentially expressed genes defined in the 4 LH subjects
    cur.idxs <- which(gene.list %in% netGenes[[choi]])
    col.idx  <- region.names == choi
    o.idx    <- region.names != choi
    
    cur.within  <- mat.dat[cur.idxs, col.idx]
    within.mean <- mean(cur.within)
    withinSE    <- sd(cur.within)/sqrt(length(cur.within))
    within.arr  <- cbind(cur.within, rep(choi, 1, length(cur.within)), rep('within', 1, length(cur.within)))
    
    if (is.null(dim(mat.dat[cur.idxs, o.idx])) == TRUE ){
      cur.betw <- mean(mat.dat[cur.idxs, o.idx])
    } else {
      cur.betw   <- rowMeans(mat.dat[cur.idxs, o.idx])
    }
    betwMean   <- mean(cur.betw)
    betwSE     <- sd(cur.betw)/sqrt(length(cur.betw))
    betw.arr   <- cbind(cur.betw, rep(choi, 1, length(cur.betw)), rep('between', 1, length(cur.betw)))
    
    aov.table <- rbind(aov.table, within.arr, betw.arr)
    dat.out   <- rbind(dat.out, c(within.mean, withinSE, betwMean, betwSE))
  }
  out <- NULL
  rownames(aov.table)  <- NULL
  aov.table            <- as.data.frame(aov.table)
  aov.table$cur.within <- as.numeric(as.character(aov.table$cur.within))
  colnames(aov.table)  <- c('log2_meanNorm_expr','network','in_out')
  dat.df <- as.data.frame(dat.out)
  rownames(dat.df)     <- region.names
  colnames(dat.df)     <- c('wmean','wse','bmean','bse')
  out[[1]] <- dat.df
  out[[2]] <- aov.table
  return(out)
}




# --------------
networkAvg <- function(regions.in, don.in, dat.type, reg2yeo, network.names, gene.arr, use.regions ){
  # Args:
  #   regions.in    = array containg indices mapping onto spatially contiguous yeo cortical parcels that will be examined
  #   don.in        = donor numbers to be examined
  #   dat.type      = field name of the data to be collapsed
  #   reg2yeo       = maps spatial indices to Yeo cortical networks
  #   network.names = the names of each Yeo cortical network  
  #   gene.arr      = list of gene names corresponding to rows in the processed expression matrix. 
  n.genes  <- length(gene.arr)
  avg.expr <- matrix(NA, n.genes, 7) # Output structure
  for (row_num in 1:n.genes){
    genes <- NULL
    for (donor in don.in){
      # Current cortical expression for an invidividual gene
      sub.expr  <- all_data[[donor]][[dat.type]][row_num, regions.in]
      cur.regs  <- reg2yeo[regions.in]
      out.means <- NULL
      # Average cortical expression within this region
      for (reg in region.names ){
        use.mean  <- mean(sub.expr[cur.regs == reg], na.rm = TRUE)
        out.means <- c(out.means, use.mean)
      }
      genes <- rbind(genes, out.means)
    }
    avg.expr[row_num,] <- colMeans(genes)
  }
  avg.expr.df <- as.data.frame(avg.expr)
  colnames(avg.expr.df) <- network.names
  rownames(avg.expr.df) <- gene.arr
  return(avg.expr.df)
}



networkAvgIndividual <- function(regions.in, don.in, dat.type, reg2yeo, network.names, gene.arr, use.regions ){
  # Args:
  #   regions.in    = array containg indices mapping onto spatially contiguous yeo cortical parcels that will be examined
  #   don.in        = donor numbers to be examined
  #   dat.type      = field name of the data to be collapsed
  #   reg2yeo       = maps spatial indices to Yeo cortical networks
  #   network.names = the names of each Yeo cortical network  
  #   gene.arr      = list of gene names corresponding to rows in the processed expression matrix. 

  dat.out <- NULL
  n.genes <- length(gene.arr)
  for (donor in don.in){
    
    avg.expr <- matrix(NA, n.genes, 7) # Output structure
    cur.data <- all_data[[donor]][[dat.type]][, regions.in]
    cur.regs <- reg2yeo[regions.in]
    ct <- 1
    for (reg in region.names ){
      use.mean      <- rowMeans(cur.data[,which(cur.regs == reg)], na.rm = TRUE)
      avg.expr[,ct] <- use.mean 
      ct <- ct + 1
    }
    colnames(avg.expr) <- region.names
    rownames(avg.expr) <- names(use.mean)
    dat.out[[donor]] <- avg.expr
  }
  
  return(dat.out)
}


replicationPlot <- function(ahba.array, ahba.sig.genes, rep.array, rep.fit, xlabel, ylabel, xlimit, ylimit){
  
  ordered.genes <- rownames(rep.fit)
  
  pos.coefs <- which(rep.fit$logFC > 0)
  q.values  <- rep.fit$adj.P.Val
  sig.coefs <- which(q.values <= .01)
  
  rep.sig.idxs  <- intersect(pos.coefs, sig.coefs)
  rep.row.names <- ordered.genes
  rep.sig.genes <- as.character(rep.row.names[rep.sig.idxs])
  
  # Color significant AHBA genes in red
  color.arr        <- rep('grey', length(ahba.array))
  match.points     <- which(toupper(ordered.genes) %in% toupper(ahba.sig.genes))
  rep.match.points <- which(toupper(ordered.genes) %in% toupper(rep.sig.genes))
  color.arr[intersect(rep.match.points, match.points)] <- 'springgreen4'
  color.arr[setdiff(rep.match.points, match.points)]   <- 'gold'
  color.arr[setdiff(match.points, rep.match.points)]   <- 'deepskyblue2'
  
  plot_df <- as.data.frame(cbind(as.numeric(ahba.array),  as.numeric(rep.array), color.arr))
  colnames(plot_df) <- c('AHBA','RepDat', 'color')
  # Keep the proper color plotting order
  plot_df$color <- factor(plot_df$color, levels = c('grey','gold','springgreen4','deepskyblue2'))
  plot_df$AHBA   <- as.numeric(as.character(plot_df$AHBA))
  plot_df$RepDat <- as.numeric(as.character(plot_df$RepDat))
  for (col in levels(plot_df$color) ){
    field.ahba   <- paste(col, 'AHBA', sep = '')
    field.repdat <- paste(col, 'RepDat', sep = '')
    
    plot_df[[field.ahba]]   <- plot_df$AHBA
    plot_df[[field.repdat]] <- plot_df$RepDat
    plot_df[[field.ahba]][which(!plot_df$color == col)] <- NA
    plot_df[[field.repdat]][which(!plot_df$color == col)] <- NA
  }
  ols <- lm(RepDat ~ AHBA,
            data = plot_df)
  ggplot( plot_df ) +
    geom_point( aes(x = greyAHBA, y = greyRepDat, color=color), size = 2.5) +
    geom_point( aes(x = goldAHBA, y = goldRepDat, color=color), size = 2.5) +
    geom_point( aes( x = deepskyblue2AHBA, y = deepskyblue2RepDat, color = color), size = 2.5) + 
    geom_point( aes(x = springgreen4AHBA, y = springgreen4RepDat, color = color), size = 2.5) +
    xlab(xlabel) + ylab(ylabel) +
    scale_color_manual(breaks=unique(plot_df$color), values=as.character(unique(plot_df$color))) + 
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
    scale_x_continuous(limits = c(-1*xlimit, xlimit), breaks = number_ticks(11))  +
    scale_y_continuous(limits = c(-1*ylimit, ylimit), breaks = number_ticks(11))  +
    geom_abline(intercept = ols$coefficients[1],
                slope = ols$coefficients[2])
}



monk.replicationPlot <- function(ahba.dat, ahba.sig.genes, rep.dat, rep.fit, xlabel, ylabel, xlimit, ylimit){
  
  ahba.array <- ahba.dat$logFC
  rep.array  <- rep.dat$logFC
  pos.coefs <- which(rep.dat$logFC > 0)
  sig.coefs <- which(rep.dat$adj.P.Val <= .01)
  rep.sig.idxs  <- intersect(pos.coefs, sig.coefs)
  rep.sig.genes <- as.character(rownames(rep.dat)[rep.sig.idxs])
  
  # Color significant AHBA genes in red
  color.arr        <- rep('grey', length(ahba.array))
  match.points     <- which(rownames(rep.dat) %in% as.character(ahba.sig.genes))
  rep.match.points <- which(rownames(rep.dat) %in% rep.sig.genes)
  
  color.arr[intersect(rep.match.points, match.points)] <- 'springgreen4'
  color.arr[setdiff(rep.match.points, match.points)]   <- 'gold'
  color.arr[setdiff(match.points, rep.match.points)]   <- 'deepskyblue2'
  
  plot_df <- as.data.frame(cbind(as.numeric(ahba.array),  as.numeric(rep.array), color.arr))
  colnames(plot_df) <- c('AHBA','RepDat', 'color')
  # Keep the proper color plotting order
  plot_df$color <- factor(plot_df$color, levels = c('grey','gold','springgreen4','deepskyblue2'))
  plot_df$AHBA   <- as.numeric(as.character(plot_df$AHBA))
  plot_df$RepDat <- as.numeric(as.character(plot_df$RepDat))
  for (col in levels(plot_df$color) ){
    field.ahba   <- paste(col, 'AHBA', sep = '')
    field.repdat <- paste(col, 'RepDat', sep = '')
    
    plot_df[[field.ahba]]   <- plot_df$AHBA
    plot_df[[field.repdat]] <- plot_df$RepDat
    plot_df[[field.ahba]][which(!plot_df$color == col)] <- NA
    plot_df[[field.repdat]][which(!plot_df$color == col)] <- NA
  }
  ols <- lm(RepDat ~ AHBA,
            data = plot_df)
  ggplot( plot_df ) +
    geom_point( aes(x = greyAHBA, y = greyRepDat, color=color), size = 2.5) +
    geom_point( aes(x = goldAHBA, y = goldRepDat, color=color), size = 2.5) +
    geom_point( aes( x = deepskyblue2AHBA, y = deepskyblue2RepDat, color = color), size = 2.5) + 
    geom_point( aes(x = springgreen4AHBA, y = springgreen4RepDat, color = color), size = 2.5) +
    xlab(xlabel) + ylab(ylabel) +
    scale_color_manual(breaks=unique(plot_df$color), values=as.character(unique(plot_df$color))) + 
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
    scale_x_continuous(limits = c(-1*xlimit, xlimit), breaks = number_ticks(11))  +
    scale_y_continuous(limits = c(-1*ylimit, ylimit), breaks = number_ticks(11))  +
    geom_abline(intercept = ols$coefficients[1],
                slope = ols$coefficients[2])
}



averageWithinCortNetworks <- function(all_data, use.regions, reg2yeo, use.donors, dat2avg, regs){
  # For each subject, get average cortical expression across all cortical regions that fall 
  #  within the same Yeo functional network
  # --------------
  # Args:
  #   all_data:    overall data struct
  #   use.regions: pre-defined array of indices for regions that meet some criteria, e.g. samples in at least 2 subs
  #   reg2yeo:     maps between region id and yeo network label
  #   use.donors:  donors you want to analyse
  #   dat2avg:     name of the field containing the data we will be averaging
  #
  # Returns:
  #     Matrix that is "# of genes" X "(# of subjects*# of regions)"
  
  ct <- 0 
  all.expr   <- NULL
  region_arr <- NULL
  donor_arr  <- NULL
  for (donor in use.donors){
    cur.cort.data <- all_data[[donor]][[dat2avg]]
    
    donor_data <- NULL
    for (reg in regs) {
      idxs_for_cur_reg   <- intersect(grep(reg, reg2yeo), use.regions)
      ct <- ct + length(idxs_for_cur_reg)
      if (is.null(dim(cur.cort.data[,idxs_for_cur_reg])) == TRUE){ # if this subject has only one region in the current network
        avgdat_for_cur_reg <- cur.cort.data[,idxs_for_cur_reg]
      } else { # more than one region
        avgdat_for_cur_reg <- rowMeans(cur.cort.data[,idxs_for_cur_reg], na.rm = TRUE)
      }
      if ( is.nan(avgdat_for_cur_reg[1]) == FALSE ){
        donor_data <- cbind(donor_data, avgdat_for_cur_reg)
        region_arr <- c(region_arr, reg)
        donor_arr  <- c(donor_arr, donor)
      }
    }
    all.expr <- cbind(all.expr, donor_data)
  }
  out <- list(all.expr, region_arr, donor_arr)
  return(out)
}


plot.blueprint <- function(blueprint.collapsed, blueprint.mature.cols, gene, regions2plot, ylimit){
  gene.idx <- which(rownames(blueprint.collapsed) %in% gene)
  gene.dat <- blueprint.collapsed[gene.idx,] - mean(blueprint.collapsed[gene.idx,])
  dat.out  <- NULL
  for (reg in regions2plot){
    reg.dat  <- gene.dat[grep(reg, names(gene.dat))]
    reg.mean <- mean(reg.dat)
    reg.sd   <- sd(reg.dat)
    reg.se   <- sd(reg.dat) / sqrt(length(reg.dat))
    dat.out  <- rbind(dat.out, c(reg.mean, reg.sd, reg.se, reg))
  }
  dat.out <- as.data.frame(dat.out)
  rownames(dat.out) <- regions2plot
  dat.out$V1 <- as.numeric(as.character(dat.out$V1))
  dat.out$V2 <- as.numeric(as.character(dat.out$V2))
  dat.out$V3 <- as.numeric(as.character(dat.out$V3))
  colnames(dat.out) <- c('mean','sd','se','reg')
  ggplot(dat.out, aes(x=reg, y=mean)) + geom_bar(stat='identity')
  dat.out$reg <- factor(dat.out$reg, levels = regions2plot)
  dat.out$mean <- as.numeric(as.character(dat.out$mean))
  # G to the G to the plot
  dodge  <- position_dodge(width=0.9)
  limits <- aes(ymax=mean+se, ymin=mean-se)
  ggplot(dat.out, aes(x=reg, y=mean)) + geom_bar(position=dodge, stat="identity") + 
    geom_errorbar(limits, position=dodge, width=0.25) +
    ylim(-ylimit, ylimit) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    ggtitle(gene)
}



averageWithinStriatNetworks <- function(all_data, donor.nums, use.striat.networks, type, net_names, atlas_field){
  # For each subject, get average striaal expression across all cortical regions that fall 
  #  within the same Choi functional network
  # --------------
  # Args:
  #   all_data:    overall data struct
  #   use.striat.networks: pre-defined array of indices for regions that meet some criteria, e.g. samples in at least 2 subs
  #   use.donors:  donors you want to analyse
  #   dat2avg:     name of the field containing the data we will be averaging
  #
  # Returns:
  #     Matrix that is "# of genes" X "(# of subjects*# of regions)"
  use.idxs    <- which(net_names %in% use.striat.networks)-1
  striat.idxs <- which(net_names %in% use.striat.networks)-1
  
  # Get Average cortical expression for each subject, in 6 
  # ------------------------------------------------------
  out.mat     <- NULL
  regions     <- NULL
  donor.array <- NULL
  for ( donor in donor.nums ) {
    cur.striat       <- all_data[[donor]][[type]]
    cur.striat.atlas <- all_data[[donor]][[atlas_field]]
    for (idx in striat.idxs){
      cdat <- cur.striat[,cur.striat.atlas == idx]
      if (is.null(dim(cdat)) == TRUE) {
        use.arr <- cdat
      } else {
        use.arr <- rowMeans(cdat)
      }
      if ( is.nan(sum(use.arr)) == FALSE ){
        out.mat     <- cbind(out.mat, use.arr)
        regions     <- c(regions, net_names[idx+1])
        donor.array <- c(donor.array, donor)
      }
    }
  }
  x <- list()
  x[[1]] <- out.mat
  x[[2]] <- regions
  x[[3]] <- donor.array
  return(x)
}



averageWithinStriatNetIndividual <- function(all_data, donor.nums, use.striat.networks, type, net_names, atlas_field){
  # For each subject, get average striaal expression across all cortical regions that fall 
  #  within the same Choi functional network
  # --------------
  # Args:
  #   all_data:    overall data struct
  #   use.striat.networks: pre-defined array of indices for regions that meet some criteria, e.g. samples in at least 2 subs
  #   use.donors:  donors you want to analyse
  #   dat2avg:     name of the field containing the data we will be averaging
  #
  # Returns:
  #     Matrix that is "# of genes" X "(# of subjects*# of regions)"
  
  use.idxs    <- which(net_names %in% use.striat.networks)-1
  striat.idxs <- which(net_names %in% use.striat.networks)-1
  
  # Get Average cortical expression for each subject, in 6 
  # ------------------------------------------------------
  all.out <- NULL
  for ( donor in donor.nums ) {
    out.mat     <- NULL
    regions     <- NULL

    cur.striat       <- all_data[[donor]][[type]]
    cur.striat.atlas <- all_data[[donor]][[atlas_field]]
    for (idx in striat.idxs){
      cdat <- cur.striat[,cur.striat.atlas == idx]
      if (is.null(dim(cdat)) == TRUE) {
        use.arr <- cdat
      } else {
        use.arr <- rowMeans(cdat)
      }
      if ( is.nan(sum(use.arr)) == FALSE ){
        out.mat     <- cbind(out.mat, use.arr)
        regions     <- c(regions, net_names[idx+1])
        colnames(out.mat) <- regions
      }
    }
    all.out[[donor]] <- out.mat
  }
  return(all.out)
}


getCortRegions <- function(all_data, filenames, atlas.num, thresh){
  all_use <- NULL
  for ( donor in filenames ){
    split.name        <- paste('splitLabel_', atlas.num, sep = '') # Each individual cortical parcel (values = 1-114)
    split.regions     <- unique(all_data[[donor]][[split.name]]) # split_region index (corresponding to the functional atlas) for each sample
    use.split.regions <- NULL
    all_use <- c(all_use, unique(split.regions[split.regions != 0]))
  }
  # regions that are present in at least X number of subjects (set with thresh variable)
  reg_counts   <- sort(table(all_use))
  regions_min  <- reg_counts[reg_counts >= thresh]
  regions.out  <- as.numeric(rownames(as.matrix(regions_min)))
  regions.out  <- sort(regions.out)
  return(regions.out)
}




query_atlas <- function(atlas, coords, rad, afni.dir) {
  cmd       <- paste(afni.dir, '/3dmaskdump -nbox ', coords[1]+rad, ':', coords[1]-rad, ' ', coords[2]+rad, ':', coords[2]-rad, ' ', coords[3]+rad, ':', coords[3]-rad, ' ', atlas, sep='')
  output    <- system(cmd, intern=TRUE)
  neigh_arr <- NULL
  for ( idx in 1:length(output) ) {
    cur_out   <- output[[idx]]
    split_out <- strsplit(cur_out, split=' ')
    network   <- as.numeric(split_out[[1]][4])
    neigh_arr <- c(neigh_arr, network)
  }
  return(neigh_arr)
}



averageCortExpr <- function(use_regions, data_struct, max_reg, sample_labels, data_type){
  # 
  # Args:
  #   use_regions: array w/ info about the network assignment of each cortical sample
  #   data.struct:    data for the current subject 
  #   max_reg:        sometimes you just want to run on the left-hemisphere, which is equivalent to atlas indices < 57
  #   sample_labels:    numeric network assignment of each sample
  #   data.type:      the name of the cortical field in 'data.struct' that we want to average
  #
  # Returns:
  #   subject-wise cortical data averaged for each network
  n_genes  <- dim(data_struct[[data_type]])[1]
  cort_arr <- matrix(NA, n_genes, max_reg)
  
  for ( reg in use_regions ){
    cur_dat <- data_struct[[data_type]][, sample_labels == reg ]
    if ( is.null(dim(cur_dat)) == TRUE ) {
      use_dat <- cur_dat
    } else {
      use_dat <- rowMeans(cur_dat)
    }
    cort_arr[,reg] <- use_dat
  }
  out           <- cort_arr
  rownames(out) <- rownames(data_struct[[data_type]])
  return(out)
}

averageStriatExpr <- function(striat.regions, data.struct, choi.names, samp.labels, data.type) {
    # 
    # Args:
    #   striat_regions: array w/ info about the network assignment of each striatal sample
    #   data.struct:    data for the current subject 
    #   choi.names:     list of all possible network assignments
    #   samp.labels:    numeric network assignment of each sample
    #   data.type:      the name of the striatal field in 'data.struct' that we want to average
    #
    # Returns:
    #   subject-wise striatal data averaged for each network
    all.expr <- NULL
    cols.out <- NULL
    for (reg in striat.regions) {
      # name array to add as column headers at the finish
      col_idx  <- samp.labels == reg # indices for this striatal subregion
      reg_data <- data.struct[[data.type]][, col_idx]
      cname    <- choi.names[as.integer(reg) + 1]
      
      # if there is only one column for this iteration
      if ( is.null(dim(reg_data)) == TRUE ){
        cur_reg <- reg_data
      } else{
        cur_reg <- rowMeans(reg_data) # Get the average across all samples
      }
      # append the data to an output array
      cols.out <- c(cols.out, cname)
      all.expr <- cbind(all.expr, cur_reg)
    }
  colnames(all.expr) <- cols.out
  striat_expr        <- as.data.frame(all.expr)
  return(striat_expr)
}


# Plot expression data from the human striatal Genotype-Phenotype Expression Project
# --------------
plot.gtex <- function(gtex.collapsed.dat.use, reg.names, gene, regions2plot, ylimit){
  gene.idx <- which(rownames(gtex.collapsed.dat.use) %in% gene)
  gene.dat <- gtex.collapsed.dat.use[gene.idx,]
  #gene.dat <- gene.dat - mean(gene.dat)
  dat.out  <- NULL
  # Calculate the mean, sd, and se for each region to be plotted
  for (reg in regions2plot){
    reg.dat  <- gene.dat[grep(reg, reg.names)]
    reg.mean <- mean(reg.dat)
    reg.sd   <- sd(reg.dat)
    reg.se   <- sd(reg.dat) / sqrt(length(reg.dat))
    dat.out  <- rbind(dat.out, c(reg.mean, reg.sd, reg.se, reg))
  }
  # Format data frame for plotting  
  dat.out <- as.data.frame(dat.out)
  rownames(dat.out) <- regions2plot
  dat.out$V1 <- as.numeric(as.character(dat.out$V1))
  dat.out$V2 <- as.numeric(as.character(dat.out$V2))
  dat.out$V3 <- as.numeric(as.character(dat.out$V3))
  colnames(dat.out) <- c('mean','sd','se','reg')
  dat.out$reg <- factor(dat.out$reg , levels=c('Nucleus_accumbens', 'Caudate', 'Putamen'))
  
  # G to the G to the plot
  dodge <- position_dodge(width=0.9)
  limits <- aes(ymax = mean + se, ymin=mean - se)
  ggplot(dat.out, aes(x=reg, y=mean)) + geom_bar(position=dodge, stat="identity") + 
    geom_errorbar(limits, position=dodge, width=0.25) +
    scale_fill_manual(values=c('dodgerblue4')) +
    ylim(-ylimit, ylimit) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    ggtitle(gene)
}


# Plot expression data from the human cortical brainspan atlas
# --------------
plot.brainspan <- function(use.eset, use_cols, gene, regions2plot, ylimit=2){
  gene.idx <- which(rownames(use.eset) %in% gene)
  gene.dat.tmp <- use.eset[gene.idx,]
  use.idxs <- which(use_cols$structure_acronym %in% regions2plot)
  gene.dat <- gene.dat.tmp[use.idxs]
  columns  <- use_cols[use.idxs,]
  dat.out  <- NULL
  names(gene.dat) <- columns$structure_acronym
  gene.dat <- gene.dat - mean(gene.dat)
  for (reg in regions2plot){
    reg.dat  <- gene.dat[grep(reg, names(gene.dat))]
    reg.mean <- mean(reg.dat)
    reg.sd   <- sd(reg.dat)
    reg.se   <- sd(reg.dat) / sqrt(length(reg.dat))
    dat.out  <- rbind(dat.out, c(reg.mean, reg.sd, reg.se, reg))
  }
  dat.out <- as.data.frame(dat.out)
  rownames(dat.out) <- regions2plot
  dat.out$V1 <- as.numeric(as.character(dat.out$V1))
  dat.out$V2 <- as.numeric(as.character(dat.out$V2))
  dat.out$V3 <- as.numeric(as.character(dat.out$V3))
  colnames(dat.out) <- c('mean','sd','se','reg')
  dat.out$reg <- factor(dat.out$reg , levels=regions2plot)
  
  dodge <- position_dodge(width=0.9)
  limits <- aes(ymax = mean + se, ymin=mean - se)
  ggplot(dat.out, aes(x=reg, y=mean)) + 
    geom_bar(position=dodge, stat="identity") + 
    geom_errorbar(limits, position=dodge, width=0.25) +
    ylim(-ylimit, ylimit) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    ggtitle(gene)
}


# Plot expression data for AHBA striatum
# --------------
ahbaPlotStriatalExpr <- function(striat.expr, striat.regions, donor.array, use.striat.networks, gene, ylimit){
  
  reg_order     <- use.striat.networks
  cur_gene_expr <- striat.expr[rownames(striat.expr) == gene]
  plot_me <- NULL
  reg_arr <- NULL
  sub_arr <- NULL
  for ( reg in reg_order ){
    reg_idxs <- grep(reg, striat.regions)
    for ( don in unique(donor.array)){
      don_idxs <- grep(don, donor.array)
      cur.val  <- cur_gene_expr[intersect(reg_idxs, don_idxs)]
      if (length(cur.val) == 0){
        cur.val <- 0
      }
      plot_me  <- c(plot_me, cur.val)
      reg_arr  <- c(reg_arr, reg)
      sub_arr  <- c(sub_arr, don)
    }
  }
  reg_fac <- factor(reg_arr, levels = reg_order)
  plot_df <- as.data.frame(cbind(plot_me, reg_arr, sub_arr))
  plot_df$reg_arr   <- factor(plot_df$reg_arr, levels = reg_order)
  colnames(plot_df) <- c('gene_expr','reg', 'sub')
  plot_df$gene_expr <- as.numeric(as.character(plot_df$gene_expr))
  
  ggplot(plot_df, aes(reg, gene_expr)) + geom_bar(aes(fill = sub), position="dodge", stat = 'identity') + 
    scale_fill_manual(values=c('gray100','gray80','gray60','gray40','gray20','gray0')) +
    scale_x_discrete( labels = levels(reg_fac)) + 
    ggtitle(gene) +
    ylim(-ylimit, ylimit) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}


# Plot expression data for AHBA cortex
# --------------
makeGenePlot <- function(all_data, gene, gene_list, n6_cort_expr, n6_cort_regions, donor_nums, ylimit){
  
  reg_order <- c('Default','Cont','Limbic','VentAttn','DorsAttn','SomMot', 'Vis')
  reg_order <- c('Limbic','Default','Cont','VentAttn','DorsAttn','SomMot', 'Vis')
  #reg_order <- c('Limbic', 'Default','Cont','VentAttn','DorsAttn','SomMot', 'Vis')
  donor_nums <- rep(donor_nums,7)
  cur_gene_expr <- n6_cort_expr[gene_list == gene]
  plot_me <- NULL
  reg_arr <- NULL
  sub_arr <- NULL
  for ( reg in reg_order ){
    reg_idxs <- grep(reg, n6_cort_regions)
    for ( don in unique(donor_nums)){
      don_idxs <- grep(don, donor_nums)
      cur.val  <- cur_gene_expr[intersect(reg_idxs, don_idxs)]
      if (length(cur.val) == 0){
        print(reg)
        print(don)
        cur.val <- 0
      }
      plot_me  <- c(plot_me, cur.val)
      reg_arr  <- c(reg_arr, reg)
      sub_arr  <- c(sub_arr, don)
    }
  }
  reg_fac <- factor(reg_arr, levels = reg_order)
  plot_df <- as.data.frame(cbind(plot_me, reg_arr, sub_arr))
  plot_df$reg_arr   <- factor(plot_df$reg_arr, levels = reg_order)
  colnames(plot_df) <- c('gene_expr','reg', 'sub')
  plot_df$gene_expr <- as.numeric(as.character(plot_df$gene_expr))
  #plot_df$gene_expr <- plot_df$gene_expr - min(plot_df$gene_expr)
  
  ggplot(plot_df, aes(reg, gene_expr)) + geom_bar(aes(fill = sub), position="dodge", stat = 'identity') + 
    scale_fill_manual(values=c('gray100','gray80','gray60','gray40','gray20','gray0')) +
    scale_x_discrete( labels = levels(reg_fac)) + 
    ggtitle(gene) + 
    ylim(-ylimit, ylimit) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}



