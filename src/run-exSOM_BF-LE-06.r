### run exSOM
### 2022-08-24_psoriasis_predicted-celltype + BF_LE_06_KS_LE_all
### 05.09.22


# load packages and source code -------------------------------------------

pacman::p_load(Seurat, tidyverse, oposSOM, here)
source(here("src","seurat-functions.r"))
source(here("src","opossom-functions.r"))
source(here("src","paths.r"))

# load data ---------------------------------------------------------------

filename <- here("2022-08-24_psoriasis_predicted-celltype.RData")
load(filename)

new_dataset_name <- "BF_LE_06_KS_LE_all"
indata_seurat_object <- load_seurat_object(new_dataset_name, output_data_path)
color_palette_seurat_object <- NULL

# run programm ------------------------------------------------------------

org_indata = env$indata
org_samples <- colnames(org_indata)
org_spot_list_group_overexpression = env$spot.list.group.overexpression
org_group_labels = env$group.labels
org_groupwise_group_colors <- env$groupwise.group.colors

# indata_seurat_object <- seurat_object
# new_dataset_name <- name_seurat_object
new_indata <- ScaleData(indata_seurat_object, rownames(indata_seurat_object)) %>% 
GetAssayData( slot = "scale.data", assay = "RNA") %>% 
as.matrix()

shared_genes <- intersect(rownames(org_indata), rownames(new_indata))

new_indata <- new_indata[shared_genes,]
new_indata <- oposSOM::Quantile.Normalization(new_indata)
new_indata <- new_indata - env$indata.gene.mean[rownames(new_indata)]

new_indata_map = matrix(NA,nrow(org_indata),ncol(new_indata),dimnames=list(rownames(org_indata),colnames(new_indata)))
new_indata_map[shared_genes,] = new_indata[shared_genes,]
# new_indata_map <- new_indata_map - mean(new_indata_map,na.rm=T)
# new_indata_map <- oposSOM::Quantile.Normalization(new_indata_map)

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

env$preferences$flip.SOM.portraits <- TRUE
if (env$preferences$flip.SOM.portraits)
{
  o <- matrix(c(1:(env$preferences$dim.1stLvlSom^2)), env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom, byrow=TRUE)
  env$som.result <- env$som.result[as.vector(o),]
}

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




