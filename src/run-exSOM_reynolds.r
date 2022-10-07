### run exSOM
### 2022-08-24_psoriasis_predicted-celltype + Reynolds2021_Data-all_seuratObj
### 07.09.22


# load packages and source code -------------------------------------------

pacman::p_load(Seurat, tidyverse, oposSOM, here)
source(here("src","seurat-functions.r"))
source(here("src","opossom-functions.r"))
source(here("src","paths.r"))

# load data ---------------------------------------------------------------

# filename <- here("oposSOM","2022-08-24_psoriasis_predicted-celltype.RData")
# load(filename)
# 
# filename <- here(output_data_path, "Reynolds2021_Data-all_seuratObj.RData")
# load(filename)
# Idents(seurat_object) <- "Cell_type_nicknames"
# indata_seurat_object <-  subset(seurat_object, downsample = 100)
# indata_seurat_object$predicted.id.Cell_type_nicknames <- indata_seurat_object$Cell_type_nicknames
# new_dataset_name <- "reynolds-2021"
# color_palette_seurat_object <- NULL

filename <- here("2022-08-24_psoriasis_predicted-celltype_exSOM_reynolds-2021 pre.RData")
load(filename)


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
    fdrtool.result <- fdrtool(env$p.g.m[,m], statistic="pvalue", verbose=FALSE, plot=FALSE)
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

3242