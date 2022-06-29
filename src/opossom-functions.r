#' Plotting gene or gene set profiles and maps of genes present in oposSOM env
#'
#' @param features character vector of features present in rownames(env$indata)
#' @param filename filename of the output pdf
#' @return None
features_profiles_and_maps <- function(features, filename) {
  filename <- paste0(filename, " - Profiles and Maps.pdf")
  pdf(filename)
  layout(matrix(c(1:8), ncol = 2, byrow = TRUE), widths = c(1, 0.5))
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
  dev.off()
}
