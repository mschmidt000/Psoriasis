#' Run scrat without seurat preprocessing and cell cycle correction
#'
#' @param ...
#' @return env

scrat_run_maria <- function(env)
{
  env$preferences$system.info <- Sys.info()
  env$preferences$session.info <- sessionInfo()
  env$preferences$started <- format(Sys.time(), "%a %d %b %Y %X")
  
  util.info("Started:", env$preferences$started)
  util.info("Name:", env$preferences$dataset.name)
  
  #### Preparation & Calculation part ####
  env <- pipeline.checkInputParameters(env)
  if (!env$passedInputChecking) {
    return(env)
  }
  
  if(env$preferences$activated.modules$reporting)
  {
    # create output directories
    dir.create(paste(env$files.name, "- Results"), showWarnings=FALSE)
    dir.create(paste(env$files.name, "- Results/CSV Sheets"), showWarnings=FALSE)
  }
  
  if(env$preferences$activated.modules$primary.analysis)
  {
    #util.info("Preprocessing Seurat Object")
    #env <- pipeline.seuratPreprocessing(env)
    
    filename <- paste(env$files.name, "pre.RData")
    util.info("Saving environment image:", filename)
    save(env, file=filename)
  }
  
  if(env$preferences$activated.modules$primary.analysis || env$preferences$activated.modules$geneset.analysis)
  {
    util.info("Loading gene annotation data.")
    env <- pipeline.prepareAnnotation(env)
  }
  
  if(env$preferences$activated.modules$primary.analysis)
  {
    #util.info("Classification of cell cycle phases")
    #env <- pipeline.cellcycleProcessing(env)
  }
  
  if (env$preferences$preprocessing$create.meta.cells)
  {
    env <- pipeline.createMetacells(env)
  }
  
  
  if(env$preferences$activated.modules$primary.analysis)
  {
    util.info("Processing SOM. This may take several time until next notification.")
    env <- pipeline.prepareIndata(env)
    env <- pipeline.generateSOM(env)
    
    filename <- paste(env$files.name, "pre.RData")
    util.info("Saving environment image:", filename)
    save(env, file=filename)
    
    util.info("Processing Differential Expression Statistics")
    env <- pipeline.calcStatistics(env)
    
    util.info("Detecting Spots")
    env <- pipeline.detectSpotsSamples(env)
    env <- pipeline.detectSpotsIntegral(env)
    env <- pipeline.patAssignment(env)
    env <- pipeline.groupAssignment(env)
  }
  
  if (env$preferences$activated.modules$geneset.analysis)
  {
    util.info("Calculating Geneset Enrichment")
    env <- pipeline.genesetStatisticSamples(env)
    env <- pipeline.genesetStatisticIntegral(env)
  }
  
  if(!is.null(env$preferences$pseudotime.estimation))
  {
    util.info("Processing Pseudotime Analysis")
    env <- pipeline.pseudotimeEstimation(env)
  }
  
  if(env$preferences$activated.modules$primary.analysis || env$preferences$activated.modules$geneset.analysis)
  {
    env$indata <- NULL
    env$indata.ensID.m <- NULL
    
    filename <- paste(env$files.name, ".RData", sep="")
    util.info("Saving environment image:", filename)
    save(env, file=filename)
    
    if (file.exists(paste(env$files.name, "pre.RData")) && file.exists(filename))
    {
      file.remove(paste(env$files.name, "pre.RData"))
    }
  }
  
  #### Reporting part ####
  
  if(env$preferences$activated.modules$reporting)
  {
    
    pipeline.summarySheetSeurat(env)
    
    if(ncol(env$metadata) < 1000)
    {
      util.info("Plotting Sample Portraits")
      pipeline.sampleExpressionPortraits(env)
    }
    
    if ( env$preferences$activated.modules$sample.similarity.analysis && ncol(env$seuratObject) > 2)
    {
      util.info("Plotting (Meta-)cell Similarity Analyses")
      dir.create(file.path(paste(env$files.name, "- Results"), "Sample Similarity Analysis"), showWarnings=FALSE)
      
      pipeline.sampleSimilarityAnalysisED(env)
      pipeline.sampleSimilarityAnalysisCor(env)
      pipeline.sampleSimilarityAnalysisICA(env)
    }
    
    util.info("Plotting Summary Sheets (Modules & PATs)")
    pipeline.summarySheetsModules(env)
    
    if(env$preferences$activated.modules$group.analysis && length(unique(env$group.labels)) >= 2)
    {
      util.info("Processing Group-centered Analyses")
      dir.create(paste(env$files.name, "- Results/Summary Sheets - Groups"), showWarnings=FALSE)
      pipeline.summarySheetsGroups(env)
    }
    
    if(!is.null(env$preferences$pseudotime.estimation))
    {
      util.info("Processing Pseudotime Reports")
      pipeline.pseudotimeReport(env)
    }
    
    util.info("Generating HTML Report")
    pipeline.htmlSummary(env)
    
  }
  
  util.info("Finished:", format(Sys.time(), "%a %b %d %X"))
  
  return(env)
}