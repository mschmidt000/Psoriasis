#' Plotting gene or gene set profiles and maps of genes present in oposSOM env
#'
#' @param env oposSOM environment
#' @param features character vector of features present in rownames(env$indata)
#' @param filename filename of the output pdf which will be saved in the results folder and extended by " - Profiles and Maps.pdf"
#' @param output_pdf outputting figure as pdf?
#' @return None
features_profiles_and_maps <- function(env, features, filename, output_pdf = TRUE) {
  
  if(output_pdf){
    filename <- here(paste(env$files.name, "- Results"), paste0(filename, " - Profiles and Maps.pdf"))
    pdf(filename)
    layout(matrix(c(1:8), ncol = 2, byrow = TRUE), widths = c(1, 0.5))
  }
  
  if (length(features) == 0) {
    return("No input features found")
  } else {
    for (i in seq_along(features)) {
      set.genes <- ifelse(is.list(features), names(features[[i]]), names(features[i]))
      if (length(intersect(set.genes, rownames(env$indata))) == 0) {
        next
      }
      feature.data <- colMeans(env$indata[intersect(set.genes, rownames(env$indata)), , drop = F], na.rm = TRUE)
      
      # barplot
      par(mar = c(1.5, 6, 1, 4))
      off.thres <- -sd(feature.data)
      on.thres <- sd(feature.data)
      barplot(feature.data,
              beside = TRUE, cex.main = 0.6,
              las = 2, cex.names = 1.2, cex.axis = 1.4, col = env$group.colors,
              border = NA,
              names.arg = rep("", ncol(env$indata))
      )
      abline(h = c(off.thres, on.thres, 0), lty = 2)
      main <- ifelse(is.list(features), "Mean log expression", "Log expression")
      mtext(main, side = 2, line = 4, cex = 0.8)
      
      # population map
      par(mar = c(2, 4, 2, 4))
      n.map <- matrix(0, env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom)
      gs.nodes <- env$som.result$nodes[set.genes]
      n.map[as.numeric(names(table(gs.nodes)))] <- table(gs.nodes)
      n.map[which(n.map == 0)] <- NA
      n.map <- matrix(n.map, env$preferences$dim.1stLvlSom)
      lim <- c(1, env$preferences$dim.1stLvlSom) + env$preferences$dim.1stLvlSom * 0.01 * c(-1, 1)
      colr <- env$color.palette.heatmaps(1000)[(na.omit(as.vector(n.map)) - min(n.map, na.rm = TRUE)) /
                                                 max(1, (max(n.map, na.rm = TRUE) - min(n.map, na.rm = TRUE))) *
                                                 999 + 1]
      plot(which(!is.na(n.map), arr.ind = TRUE),
           xlim = lim, ylim = lim, pch = 16, axes = FALSE,
           xlab = "", ylab = "", xaxs = "i", yaxs = "i", col = colr,
           cex = 0.5 + na.omit(as.vector(n.map)) / max(n.map, na.rm = TRUE) * 1.8
      )
      main <- ifelse(is.list(features), names(features)[i], features[i])
      mtext(main, side = 2, line = 2.5, cex = 0.5)
      box()
    }
  }
  if(output_pdf){
    dev.off()
  }
}

#' Resort oposSOM env by new labels
#'
#' @param env oposSOM enviroment
#' @param new_labels character vector with new labels and sorted by colnames(env$indata)
#' @param group_colors color palette named by group.labels
#' @return env
resort_env <- function(env = env, new_labels = new_labels, labels_name = labels_name, group_colors = NULL){
  
  env$preferences$activated.modules$primary.analysis <- FALSE
  env$preferences$activated.modules$geneset.analysis <- TRUE
  env$files.name <- env$preferences$dataset.name <- paste(env$preferences$dataset.name, labels_name, sep = "_")
  
  env$preferences$system.info <- Sys.info()
  env$preferences$session.info <- sessionInfo()
  env$preferences$started <- format(Sys.time(), "%a %d %b %Y %X")
  
  new_order <- names(sort(new_labels))
  
  env$group.labels <- new_labels[new_order]
  env$indata <- env$indata[,new_order]
  if(is.null(group_colors)){
  env$groupwise.group.colors <- rainbow(length(unique(env$group.labels)))
  names(env$groupwise.group.colors) <- unique(env$group.labels)
  } else {
    env$groupwise.group.colors <- group_colors
    env$groupwise.group.colors <- env$groupwise.group.colors[unique(env$group.labels)]
  }
  env$group.colors <- env$groupwise.group.colors[match(env$group.labels, names(env$groupwise.group.colors))]
  
  return(env)
}

#' Rerun necessary oposSOM analysis-pipeline parts and output report
#'
#' @param env oposSOM enviroment
#' @return NULL
analysis_and_output_after_resort <- function(env, only_output = FALSE){
  
  if(!only_output){
  
    util.info("Processing Differential Expression Statistics")
    env <- pipeline.diffExpressionStatistics(env)
    
    util.info("Detecting Spots")
    env <- pipeline.detectSpotsSamples(env)
    env <- pipeline.detectSpotsModules(env)
    env <- pipeline.patAssignment(env)
    env <- pipeline.groupAssignment(env) 
    
    if (env$preferences$activated.modules$geneset.analysis)
    {
      util.info("Calculating Geneset Enrichment")
      env <- pipeline.genesetStatisticSamples(env)
      env <- pipeline.genesetStatisticModules(env)
    }
    
    if (env$preferences$activated.modules$psf.analysis)
    {
      util.info("Calculating Pathway Signal Flow (PSF)")
      env <- pipeline.PSFcalculation(env)    
    }
    
    if(env$preferences$activated.modules$primary.analysis || env$preferences$activated.modules$geneset.analysis)
    {    
      filename <- paste(env$files.name, ".RData", sep="")
      util.info("Saving environment image:", filename)
      save(env, file=filename)
      
      if (file.exists(paste(env$files.name, "pre.RData")) && file.exists(filename))
      {
        file.remove(paste(env$files.name, "pre.RData"))
      }
    }  
  
  }
  
  #### Reporting part ####
  
  if(env$preferences$activated.modules$reporting)
  {
    
    util.info("Plotting Supporting Information")
    pipeline.supportingMaps(env)
    pipeline.entropyProfiles(env)
    pipeline.topologyProfiles(env)
    
    
    if(length(env$chromosome.list) > 0)
    {
      util.info("Plotting Chromosome Expression Reports")
      pipeline.chromosomeExpressionReports(env)
    }
    
    if(ncol(env$indata) < 1000)
    {
      util.info("Plotting Sample Portraits")
      pipeline.sampleExpressionPortraits(env)
    } 
    
    if ( env$preferences$activated.modules$sample.similarity.analysis && ncol(env$indata) > 2)
    {    
      util.info("Plotting Sample Similarity Analysis")
      dir.create("Sample Similarity Analysis", showWarnings=FALSE)
      
      pipeline.sampleSimilarityAnalysisED(env)
      pipeline.sampleSimilarityAnalysisCor(env)
      pipeline.sampleSimilarityAnalysisICA(env)
      pipeline.sampleSimilarityAnalysisSOM(env)
    }
    
    if (env$preferences$activated.modules$geneset.analysis)
    {
      dir.create("Geneset Analysis", showWarnings=FALSE)
      
      util.info("Plotting Geneset Enrichment Heatmaps")
      pipeline.genesetOverviews(env)
      
      util.info("Plotting Geneset Profiles and Maps")
      pipeline.genesetProfilesAndMaps(env)
      
      util.info("Calculating Cancer Hallmark Enrichment")
      pipeline.cancerHallmarks(env)
    }
    
    if (env$preferences$activated.modules$psf.analysis)
    {
      util.info("Plotting PSF results")
      pipeline.PSFoutput(env)
    }
    
    util.info("Writing Gene Lists")
    pipeline.geneLists(env)
    
    util.info("Plotting Summary Sheets (Samples)")
    pipeline.summarySheetsSamples(env)
    
    util.info("Plotting Summary Sheets (Modules & PATs)")
    pipeline.summarySheetsModules(env)
    pipeline.summarySheetsPATs(env)
    
    
    if(env$preferences$activated.modules$group.analysis && length(unique(env$group.labels)) >= 2)
    {
      util.info("Processing Group-centered Analyses")
      pipeline.groupAnalysis(env)
    }
    
    if(env$preferences$activated.modules$difference.analysis)
    {
      util.info("Processing Difference Analyses")
      pipeline.differenceAnalyses(env)
    }
    
    util.info("Generating HTML Report")
    pipeline.htmlSampleSummary(env)
    pipeline.htmlModuleSummary(env)
    pipeline.htmlGenesetAnalysis(env)  
    pipeline.htmlPsfAnalysis(env)
    pipeline.htmlSummary(env)
    
  }    
  
  util.info("Finished:", format(Sys.time(), "%a %b %d %X"))
}

#' Creating opossom env for samples being groups
#'
#' @param env oposSOM enviroment
#' @return local.env
create_local_env_for_groups <- function(env){
  
  local.env <- new.env()
  local.env$preferences <- env$preferences
  local.env$gene.info <- env$gene.info
  local.env$gs.def.list <- env$gs.def.list
  local.env$som.result <- env$som.result
  local.env$files.name <- env$files.name
  local.env$csv.function <- env$csv.function
  local.env$color.palette.portraits <- env$color.palette.portraits
  local.env$color.palette.heatmaps <- env$color.palette.heatmaps
  
  # calculate differential expression statistics
  
  local.env$p.g.m <- matrix(NA, nrow(env$indata), length(unique(env$group.labels)),
                            dimnames=list(rownames(env$indata), unique(env$group.labels)))
  local.env$fdr.g.m <- matrix(NA, nrow(env$indata), length(unique(env$group.labels)),
                              dimnames=list(rownames(env$indata), unique(env$group.labels)))
  local.env$n.0.m <- rep(NA, length(unique(env$group.labels)))
  names(env$n.0.m) <- unique(env$group.labels)
  local.env$perc.DE.m <- rep(NA, length(unique(env$group.labels)))
  names(env$perc.DE.m) <- unique(env$group.labels)
  
  
  for (gr in seq_along(unique(env$group.labels)))
  {
    samples.indata <- which(env$group.labels==unique(env$group.labels)[gr])
    
    local.env$p.g.m[,gr] <- apply( env$indata, 1, function(x)
    {
      if( var(x[-samples.indata]) == 0 ) return(1) 
      p <- t.test( x[samples.indata], x[-samples.indata], var.equal=length(samples.indata)==1 )$p.value
      if( p < 1e-16) p <- 1e-16
      return( p )
    } )
    
    suppressWarnings({
      try.res <- try({
        fdrtool.result <- fdrtool(local.env$p.g.m[,gr], statistic="pvalue", verbose=FALSE, plot=FALSE)
      }, silent=TRUE)
    })
    
    if (!is(try.res,"try-error"))
    {
      #      p.g.m[,gr] <- fdrtool.result$pval
      local.env$fdr.g.m[,gr] <- fdrtool.result$lfdr
      local.env$n.0.m[gr] <- fdrtool.result$param[1,"eta0"]
      local.env$perc.DE.m[gr] <- 1 - local.env$n.0.m[gr]
    } else
    {
      #      p.g.m[,gr] <- order(apply(indata[,samples.indata,drop=FALSE],1,mean)) / nrow(indata)
      local.env$fdr.g.m[,gr] <- local.env$p.g.m[,gr]
      local.env$n.0.m[gr] <- 0.5
      local.env$perc.DE.m[gr] <- 0.5
    }    
  }
  
  
  
  
  # average over group members
  
  local.env$metadata <- do.call(cbind, by(t(env$metadata), env$group.labels, colMeans)[unique(env$group.labels)])
  
  local.env$indata <- do.call(cbind, by(t(env$indata+env$indata.gene.mean),
                                        env$group.labels,
                                        colMeans)[unique(env$group.labels)])
  
  local.env$indata.gene.mean <- rowMeans(local.env$indata)
  
  if (local.env$preferences$feature.centralization)
  {
    local.env$indata <- local.env$indata - local.env$indata.gene.mean
  }
  
  local.env$group.colors <- env$group.colors[match(colnames(local.env$indata), env$group.labels)]
  local.env$group.labels <- env$group.labels[match(colnames(local.env$indata), env$group.labels)]
  names(local.env$group.labels) <- local.env$group.labels
  names(local.env$group.colors) <- local.env$group.labels
  
  
  
  local.env$output.paths <- c("CSV" = "Summary Sheets - Groups/CSV Sheets",
                              "Summary Sheets Samples"= "Summary Sheets - Groups/Reports")
  
  local.env <- pipeline.detectSpotsSamples(local.env)
  
  if (local.env$preferences$activated.modules$geneset.analysis)
  {
    local.env <- pipeline.genesetStatisticSamples_adapted(local.env)
  }
  
  return(local.env)
}

#' Transfer a color palette to pastel colors
#'
#' @param palette color palette to transfer to pastel colors
#' @param lightness lightness of colors, number between 0 and 1
#' @return character vector with pastel colors
make_pastel_colors <- function(palette, lightness){
  a <- palette
  a1 <- col2rgb(a)
  a2 <- rgb2hsv(a1)
  n <- lightness
  return(hsv(a2[1,], a2[2,]*n, a2[3,]))
}

#' Run extension SOM based on oposSOM
#'
#' @param env oposSOM environment to extend
#' @param seurat_object data set by which the oposSOM env should be extended
#' @return name of the new data set by which the oposSOM env should be extended
run_exSOM <- function(env, seurat_object, name_seurat_object, color_palette_seurat_object = NULL){
  
  org_indata = env$indata
  org_samples <- colnames(org_indata)
  org_spot_list_group_overexpression = env$spot.list.group.overexpression
  org_group_labels = env$group.labels
  org_groupwise_group_colors <- env$groupwise.group.colors
  
  indata_seurat_object <- seurat_object
  new_dataset_name <- name_seurat_object
  
  new_indata <- as.matrix(GetAssayData(indata_seurat_object, slot = "scale.data", assay = "RNA"))
  
  shared_genes <- intersect(rownames(org_indata), rownames(new_indata))
  
  new_indata <- new_indata[shared_genes,]
  new_indata <- oposSOM::Quantile.Normalization(new_indata)
  new_indata <- new_indata - env$indata.gene.mean[rownames(new_indata)]
  
  new_indata_map = matrix(NA,nrow(org_indata),ncol(new_indata),dimnames=list(rownames(org_indata),colnames(new_indata)))
  new_indata_map[shared_genes,] = new_indata[shared_genes,]
  new_indata_map <- new_indata_map - mean(new_indata_map,na.rm=T)
  new_indata_map <- oposSOM::Quantile.Normalization(new_indata_map)
  
  env$files.name <- paste(env$files.name, "exSOM", new_dataset_name, sep = "_")
  env$preferences$dataset.name <- env$files.name
  env$preferences$system.info <- Sys.info()
  env$preferences$session.info <- sessionInfo()
  env$preferences$started <- format(Sys.time(), "%a %d %b %Y %X")
  
  new_labels <- as.character(indata_seurat_object$predicted.id.Cell_type_nicknames)
  new_labels <- paste(new_labels, "ex")
  names(new_labels) <- rownames(indata_seurat_object@meta.data)
  new_samples_sort <- names(sort(new_labels))
  new_labels <- new_labels[new_samples_sort]
  
  indata_combine = cbind( org_indata, new_indata_map[rownames(org_indata),new_samples_sort])
  
  env$group.labels = c(org_group_labels, new_labels)
  
  env$indata = indata_combine
  
  new_groupwise_labels <- unique(new_labels)
  
  if(is.null(color_palette_seurat_object)){
    
    new_groupwise_colors <- make_pastel_colors(rainbow(length(new_groupwise_labels)), 0.6)
    names(new_groupwise_colors) <- new_groupwise_labels
    
    env$groupwise.group.colors <- c(org_groupwise_group_colors, new_groupwise_colors)
    
  } else {
    
    names(color_palette_seurat_object) <- paste(names(color_palette_seurat_object), "ex")
    color_palette_seurat_object <- color_palette_seurat_object[names(new_groupwise_labels)]
    env$groupwise.group.colors <- c(org_groupwise_group_colors, color_palette_seurat_object)
    
  }
  env$group.colors = env$groupwise.group.colors[match(env$group.labels, names(env$groupwise.group.colors))]
  names(env$group.colors) <- names(env$group.labels)
  
  batch <- rep("1",ncol(env$indata))
  batch[which(colnames(env$indata)  %in% colnames(new_indata))] <- "2"
  env$som.result <- oposSOM::som.linear.init(indata_combine,somSize=env$preferences$dim.1stLvlSom, batch=batch)
  env$som.result[,which(colnames(env$indata)  %in% colnames(new_indata))]=0

  env$som.result <- oposSOM::som.training( indata_combine, env$som.result, metricSamples= which(colnames(env$som.result)  %in% colnames(org_indata)), verbose = TRUE )
  env$metadata <- env$som.result$weightMatrix
  colnames(env$metadata) <- colnames(env$indata)
  env$som.result$weightMatrix <- NULL
  
  dir.create(paste(env$files.name, "- Results"), showWarnings=FALSE)
  dir.create(paste(env$files.name, "- Results/CSV Sheets"), showWarnings=FALSE)
  
  filename <- paste(env$files.name, "pre.RData")
  util.info("Saving environment image:", filename)
  save(env, file=filename)
  
  var <- stats::var
  env <- pipeline.diffExpressionStatistics(env)
  env <- pipeline.detectSpotsSamples(env)
  env <- pipeline.detectSpotsModules(env)
  
  filename <- paste(env$files.name, "pre.RData")
  util.info("Saving environment image:", filename)
  save(env, file=filename)
  
  
  env$spot.list.group.overexpression$overview.map =  org_spot_list_group_overexpression$overview.map
  env$spot.list.group.overexpression$overview.mask =  org_spot_list_group_overexpression$overview.mask
  env$spot.list.group.overexpression$spots =  org_spot_list_group_overexpression$spots
  env$spot.list.group.overexpression$spotdata <- t(sapply(env$spot.list.group.overexpression$spots, function(x)
  {
    if (length(x$genes > 0))
    {
      colMeans(env$indata[x$genes,,drop=FALSE],na.rm=TRUE)
    } else
    {
      rep(0, ncol(env$indata))
    }
  }))
  colnames(env$spot.list.group.overexpression$spotdata) <- colnames(env$indata)
  
  env <- pipeline.patAssignment(env)
  env <- pipeline.groupAssignment(env)
  env$preferences$activated.modules$geneset.analysis=F
  
  ### Reporting part ###
  
  filename <- paste(env$files.name, ".RData", sep="")
  util.info("Saving environment image:", filename)
  save(env, file=filename)
  
  if (file.exists(paste(env$files.name, "pre.RData")) && file.exists(filename))
  {
    file.remove(paste(env$files.name, "pre.RData"))
  }

  util.info("Plotting Supporting Information")
  dir.create(paste(env$files.name, "- Results/Supporting Maps&Profiles"), showWarnings=FALSE)
  pipeline.supportingMaps(env)
  pipeline.entropyProfiles(env)
  pipeline.topologyProfiles(env)

  if(length(env$chromosome.list) > 0)
  {
    util.info("Plotting Chromosome Expression Reports")
    pipeline.chromosomeExpressionReports(env)
  }
  
  if(ncol(env$indata) < 1000)
  {
    util.info("Plotting Sample Portraits")
    pipeline.sampleExpressionPortraits(env)
  } 
  
  if ( env$preferences$activated.modules$sample.similarity.analysis && ncol(env$indata) > 2)
  {    
    util.info("Plotting Sample Similarity Analysis")
    dir.create("Sample Similarity Analysis", showWarnings=FALSE)
    
    pipeline.sampleSimilarityAnalysisED(env)
    pipeline.sampleSimilarityAnalysisCor(env)
    pipeline.sampleSimilarityAnalysisICA(env)
    pipeline.sampleSimilarityAnalysisSOM(env)
  }
  
  if (env$preferences$activated.modules$geneset.analysis)
  {
    dir.create("Geneset Analysis", showWarnings=FALSE)
    
    util.info("Plotting Geneset Enrichment Heatmaps")
    pipeline.genesetOverviews(env)
    
    util.info("Plotting Geneset Profiles and Maps")
    pipeline.genesetProfilesAndMaps(env)
    
    util.info("Calculating Cancer Hallmark Enrichment")
    pipeline.cancerHallmarks(env)
  }
  
  if (env$preferences$activated.modules$psf.analysis)
  {
    util.info("Plotting PSF results")
    pipeline.PSFoutput(env)
  }
  
  util.info("Writing Gene Lists")
  pipeline.geneLists(env)
  
  util.info("Plotting Summary Sheets (Samples)")
  pipeline.summarySheetsSamples(env)
  
  util.info("Plotting Summary Sheets (Modules & PATs)")
  pipeline.summarySheetsModules(env)
  pipeline.summarySheetsPATs(env)
  
  
  if(env$preferences$activated.modules$group.analysis && length(unique(env$group.labels)) >= 2)
  {
    util.info("Processing Group-centered Analyses")
    pipeline.groupAnalysis(env)
  }
  
  if(env$preferences$activated.modules$difference.analysis)
  {
    util.info("Processing Difference Analyses")
    pipeline.differenceAnalyses(env)
  }
  
  util.info("Generating HTML Report")
  pipeline.htmlSampleSummary(env)
  pipeline.htmlModuleSummary(env)
  pipeline.htmlGenesetAnalysis(env)  
  pipeline.htmlPsfAnalysis(env)
  pipeline.htmlSummary(env)
  

util.info("Finished:", format(Sys.time(), "%a %b %d %X"))
  
filename <- paste(env$files.name, ".RData", sep="")
save(env, file=filename)

return(env)
}

#' Creating opossom env for samples being differences between specific groups
#'
#' @param env oposSOM environment
#' @param differences.list list with group.labels for which differential expression should be calculated
#' @return env
create_local_env_for_differences <- function(env, differences.list){
  
  local.env <- new.env()
  local.env$preferences <- env$preferences
  local.env$gene.info <- env$gene.info
  local.env$gs.def.list <- env$gs.def.list
  local.env$som.result <- env$som.result
  local.env$files.name <- env$files.name
  local.env$csv.function <- env$csv.function
  local.env$color.palette.portraits <- env$color.palette.portraits
  local.env$color.palette.heatmaps <- env$color.palette.heatmaps
  local.env$indata.gene.mean <- env$indata.gene.mean
  
  local.env$p.g.m <- matrix(NA, nrow(env$indata), length(differences.list),
                          dimnames=list(rownames(env$indata), names(differences.list)))
  
  local.env$fdr.g.m <- matrix(NA, nrow(env$indata), length(differences.list),
                            dimnames=list(rownames(env$indata), names(differences.list)))
  
  local.env$n.0.m <- rep(NA, length(differences.list))
  names(local.env$n.0.m) <- names(differences.list)
  
  local.env$perc.DE.m <- rep(NA, length(differences.list))
  names(local.env$perc.DE.m) <- names(differences.list)
  
  indata.d <- matrix(NA, nrow(env$indata), length(differences.list),
                   dimnames=list(rownames(env$indata), names(differences.list)))
  metadata.d <- matrix(NA, nrow(env$metadata), length(differences.list),
                     dimnames=list(rownames(env$metadata), names(differences.list)))
  
  
  env2 <- env
  load("/scratch/maria/Wound_Healing_Ruben_Ferrer/Haensel_2019/UW&PBS&Sit_oposSOMinput_group.labels.colors_500.RData")
  for(l in unique(group.labels.colors[[1]])){
  
  for (d in seq_along(differences.list)){
    env <- env2
    samples.cluster <- names(group.labels.colors[[1]])[which(group.labels.colors[[1]] == l)]
    samples.indata <-
      list(intersect(names(differences.list[[d]][[2]]), samples.cluster), intersect(names(differences.list[[d]][[1]]), samples.cluster))
    
    
    indata.d[,d] <- rowMeans(env$indata[,samples.indata[[1]],drop=FALSE]) -
      rowMeans(env$indata[,samples.indata[[2]],drop=FALSE])
    
    metadata.d[,d] <- rowMeans(env$metadata[,samples.indata[[1]],drop=FALSE]) -
      rowMeans(env$metadata[,samples.indata[[2]],drop=FALSE])
    
    local.env$p.g.m[,d] <- apply( env$indata, 1, function(x)
    {
      if( length(samples.indata[[1]])>1 && var(x[samples.indata[[1]]]) == 0 ) return(1) 
      if( length(samples.indata[[2]])>1 && var(x[samples.indata[[2]]]) == 0 ) return(1) 
      
      return( t.test( x[samples.indata[[1]]], x[samples.indata[[2]]], paired=FALSE, var.equal=(length(samples.indata[[1]])==1 || length(samples.indata[[2]])==1 ) )$p.value )
    } )
    
    suppressWarnings({
      try.res <- try({
        fdrtool.result <- fdrtool(local.env$p.g.m[,d], statistic="pvalue", verbose=FALSE, plot=FALSE)
      }, silent=TRUE)
    })
    
    if (!is(try.res,"try-error"))
    {
      local.env$fdr.g.m[,d] <- fdrtool.result$lfdr
      local.env$n.0.m[d] <- fdrtool.result$param[1,"eta0"]
      local.env$perc.DE.m[d] <- 1 - local.env$n.0.m[d]
    } else
    {
      local.env$fdr.g.m[,d] <- local.env$p.g.m[,d]
      local.env$n.0.m[d] <- 0.5
      local.env$perc.DE.m[d] <- 0.5
    }
    
  }
  
  local.env$indata <- indata.d
  colnames(local.env$indata) <- names(differences.list)
  
  local.env$metadata <- metadata.d
  colnames(local.env$metadata) <- names(differences.list)
  
  local.env$group.labels <- names(differences.list)
  names(local.env$group.labels) <- names(differences.list)
  local.env$group.colors <- rep("gray20",length(differences.list))
  names(local.env$group.colors) <- names(differences.list)
  
  local.env$output.paths <- c("CSV" = paste(env$files.name, "- Results/Summary Sheets - Differences/CSV Sheets"),
                              "Summary Sheets Samples"= paste(env$files.name, "- Results/Summary Sheets - Differences/Reports"))
  
  local.env <- pipeline.detectSpotsSamples(local.env)
  
  if (local.env$preferences$activated.modules$geneset.analysis)
  {
    if (ncol(local.env$indata) == 1)   # crack for by command, which requires >=2 columns
    {
      local.env$indata <- cbind(local.env$indata, local.env$indata)
      local.env <- pipeline.detectSpotsSamples(local.env)
      # local.env <- pipeline.genesetStatisticSamples(local.env)  
      ### perform GS analysis ###
      local.env$indata.ensID.m <- local.env$indata[local.env$gene.info$ensembl.mapping[,1],]
      local.env$indata.ensID.m <- do.call(rbind, by(local.env$indata.ensID.m, local.env$gene.info$ensembl.mapping[,2], colMeans))
      
      mean.ex.all <- colMeans( local.env$indata.ensID.m )
      sd.ex.all <- apply( local.env$indata.ensID.m, 2, sd )
      
      gs.null.list <- list()
      for (i in seq_along(local.env$gs.def.list))
      {
        gs.null.list[[i]] <-
          list(Genes=sample(unique(local.env$gene.info$ensembl.mapping$ensembl_gene_id), length(local.env$gs.def.list[[i]]$Genes)))
      }
      
      null.scores <- sapply( gs.null.list, Sample.GSZ, local.env$indata.ensID.m, mean.ex.all, sd.ex.all )
      null.culdensity <- ecdf(abs(unlist(null.scores)))
      
      local.env$samples.GSZ.scores <- t( sapply( local.env$gs.def.list, Sample.GSZ, local.env$indata.ensID.m, mean.ex.all, sd.ex.all ) )
      
      local.env$spot.list.samples <- lapply(seq_along(local.env$spot.list.samples) , function(m)
      {
        x <- local.env$spot.list.samples[[m]]
        
        x$GSZ.score <- local.env$samples.GSZ.scores[,m]
        x$GSZ.p.value <- 1 - null.culdensity(abs(x$GSZ.score))
        names(x$GSZ.p.value) <- names(x$GSZ.score)
        
        return(x)
      })
      names(local.env$spot.list.samples) <- colnames(local.env$indata)
      local.env$indata <- local.env$indata[,1,drop=FALSE]
      local.env$spot.list.samples <- local.env$spot.list.samples[1]
      local.env$samples.GSZ.scores <- local.env$samples.GSZ.scores[,1,drop=FALSE]
    } else
    {
      local.env <- pipeline.genesetStatisticSamples(local.env)
    }
    }
  }
  return(env)
}



pipeline.diffExpressionStatistics_adapted <- function(env)
{
  util.info("Calculating Single Gene Statistic")
  progressbar <- newProgressBar(min = 0, max = ncol(env$indata)); cat("\r")
  
  ### Single Genes ###
  
  env$p.g.m <- matrix(NA, nrow(env$indata), ncol(env$indata), dimnames=list(rownames(env$indata), colnames(env$indata)))
  
  env$n.0.m <- rep(NA, ncol(env$indata))
  names(env$n.0.m) <- colnames(env$indata)
  
  env$perc.DE.m <- rep(NA, ncol(env$indata))
  names(env$perc.DE.m) <- colnames(env$indata)
  
  env$fdr.g.m <- matrix(NA, nrow(env$indata), ncol(env$indata), dimnames=list(rownames(env$indata), colnames(env$indata)))
  
  
  for (m in 1:ncol(env$indata))
  {
    env$p.g.m[,m] <- apply( env$indata, 1, function(x)
    {
      if( all(x[-m] == x[-m][1]) ) return(1) 
      
      return( t.test( x[m], x[-m], var.equal=TRUE )$p.value )
    } )
    
    suppressWarnings({ try.res <- try({
      fdrtool.result <- fdrtool(na.omit( local.env$indata.ensID.m), statistic="pvalue", verbose=FALSE, plot=FALSE)
    }, silent=TRUE) })
    
    if (!is(try.res,"try-error"))
    {
      env$fdr.g.m[,m] <- fdrtool.result$lfdr
      
      env$n.0.m[m] <- fdrtool.result$param[1,"eta0"]
      env$perc.DE.m[m] <- 1 - env$n.0.m[m]
    } else # happens for eg phenotype data
    {
      env$p.g.m[,m] <- order(env$indata[,m]) / nrow(env$indata)
      env$fdr.g.m[,m] <- env$p.g.m[,m]
      
      env$n.0.m[m] <- 1
      env$perc.DE.m[m] <- 0
    }
    
    setTxtProgressBar( progressbar, progressbar$getVal()+1 )
  }
  progressbar$kill()
  
  
  ### Metagenes ###
  
  util.info("Calculating Metagene Statistic")
  
  env$p.m <- matrix(NA, env$preferences$dim.1stLvlSom ^ 2, ncol(env$indata),
                    dimnames=list(1:(env$preferences$dim.1stLvlSom ^ 2), colnames(env$indata)))
  
  p.m.help <- do.call(rbind, by(env$p.g.m, env$som.result$feature.BMU, colMeans))
  env$p.m[rownames(p.m.help),] <- p.m.help
  
  return(env)
}


#' Creating opossom env for samples being differences between specific groups
#'
#' @param env oposSOM environment
#' @param differences.list list with group.labels for which differential expression should be calculated
#' @return env
flip_SOM <- function(env, differences.list){
  temp_matrix <- matrix(1:env$preferences$dim.1stLvlSom^2, ncol = env$preferences$dim.1stLvlSom, byrow = TRUE)
  temp_matrix_rotated <- temp_matrix[env$preferences$dim.1stLvlSom:1,]
  temp_vector_rotated <- as.vector(t(temp_matrix_rotated))
  names(temp_vector_rotated) <- 1:env$preferences$dim.1stLvlSom^2
  
  
  env$metadata <- env$metadata[env$preferences$dim.1stLvlSom:1,]
  
  env$som.result$node.summary$y <- rep(env$preferences$dim.1stLvlSom:1, each = env$preferences$dim.1stLvlSom)
  env$som.result$node.summary <- env$som.result$node.summary[order(env$som.result$node.summary$y),]
  
  for(i in 1:env$preferences$dim.1stLvlSom^2){
    if(length(length(env$som.result$feature.BMU[which(env$som.result$feature.BMU %in% as.numeric(names(temp_vector_rotated))[i])])>0)){
      env$som.result$feature.BMU[which(env$som.result$feature.BMU %in% as.numeric(names(temp_vector_rotated))[i])] <- temp_vector_rotated[i]
    }
  }

  env$gene.info$coordinates <- apply( env$som.result$node.summary[env$som.result$feature.BMU,c("x","y")], 1, paste, collapse=" x " )
  names(env$gene.info$coordinates) <- rownames(env$indata)
  
  
}