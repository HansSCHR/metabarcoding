amp_rarecurve <- function(data, step = 100, ylim = NULL, xlim = NULL, label = F, color = NULL, legend = T, color.vector = NULL, legend.position = "topleft"){
  
  abund = otu_table(data)@.Data %>% as.data.frame()
  
  if (!is.null(color)) {
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length=n+1)
      hcl(h=hues, l=65, c=100)[1:n]
    }
    group_vector<-sample_data(data)[,color]@.Data %>% as.data.frame()
    names(group_vector)<-"color_variable"
    group_vector<-as.character(group_vector$color_variable)
    groups<-unique(group_vector)
    n = length(groups)
    cols = gg_color_hue(n)
    if (!is.null(color.vector)){ cols <- color.vector}
    
    col_vector<-rep("black",length(group_vector))
    for (i in 1:length(group_vector)){
      col_vector[i]<-cols[match(group_vector[i],groups)]
    }
  } else {
    col_vector = "black"
  }
  
  if (is.null(ylim) & is.null(xlim)){
    rarecurve(t(abund), step = step, label = label, col = col_vector)
  }
  if (!is.null(ylim) & !is.null(xlim)){
    rarecurve(t(abund), step = step, ylim = ylim, xlim = xlim, label = label, col = col_vector)
  }
  if (!is.null(ylim) & is.null(xlim)){
    rarecurve(t(abund), step = step, ylim = ylim, label = label, col = col_vector)
  }
  if (is.null(ylim) & !is.null(xlim)){
    rarecurve(t(abund), step = step, xlim = xlim, label = label, col = col_vector)
  }
  
  if (!is.null(color) & legend == T){
    legend(legend.position,legend = groups,fill = cols, bty = "n")
  }
}

plot_composition <- function(physeq, taxaRank1 = "Phylum", taxaSet1 = "Proteobacteria",
                             taxaRank2 = "Family", numberOfTaxa = 9, fill = NULL,
                             x = "Sample", y = "Abundance", facet_grid = NULL) {
  ## Args:
  ## - physeq: phyloseq class object
  ## - taxaRank1: taxonomic level in which to do the first subsetting
  ## - taxaSet1: subset of level taxaRank1 to use
  ## - taxaRank2: taxonomic level used to agglomerate
  ## - numberOfTaxa: number of (top) taxa to keep at level taxaRank2
  ##
  ## Returns:
  ## - ggplot2 graphics
  ggdata <- ggformat(physeq, taxaRank1, taxaSet1, taxaRank2, numberOfTaxa)
  p <- ggplot(ggdata, aes_string(x = x, y = y, fill = fill, color = fill, group = "Sample"))
  ## Manually change color scale to assign grey to "Unknown" (if any)
  if (!is.null(fill) && any(c("Unknown", "Other") %in% unique(ggdata[, fill]))) {
    ranks <- as.character(unique(ggdata[, fill]))
    ranks <- ranks[ ! ranks %in% c("Unknown", "Other")]
    colvals <- c(gg_color_hue(length(ranks)), "grey45", "black")
    names(colvals) <- c(ranks, "Unknown", "Other")
    ## Now add the manually re-scaled layer with Unassigned as grey
    p <- p + scale_fill_manual(values=colvals) + scale_color_manual(values = colvals)
    
  }
  p <- p + geom_bar(stat = "identity", position = "stack")
  if ( !is.null(facet_grid)) {
    p <- p + facet_grid(facets = facet_grid, scales = "free_x")
  }
  p <- p + theme(axis.text.x=element_text(angle=90), axis.title.x=element_blank())
  p <- p + ggtitle(paste("Composition within", taxaSet1, "(", numberOfTaxa, "top", taxaRank2, ")"))
  return(p)
}

ggrare <- function(physeq_object, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  
  x <- methods::as(phyloseq::otu_table(physeq_object), "matrix")
  if (phyloseq::taxa_are_rows(physeq_object)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  # Get sample data
  if (!is.null(phyloseq::sample_data(physeq_object, FALSE))) {
    sdf <- methods::as(phyloseq::sample_data(physeq_object), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  # Add, any custom-supplied plot-mapped variables
  if ( length(color) > 1 ) {
    data$color <- color
    names(data)[names(data) == "color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  
  if ( length(label) > 1 ) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot2::ggplot(data = data,
                       ggplot2::aes_string(x = "Size",
                                           y = ".S",
                                           group = "Sample",
                                           color = color))
  
  p <- p + ggplot2::labs(x = "Sequence Sample Size", y = "Species Richness")
  
  if (!is.null(label)) {
    p <- p + ggplot2::geom_text(data = labels,
                                ggplot2::aes_string(x = "x",
                                                    y = "y",
                                                    label = label,
                                                    color = color),
                                size = 4, hjust = 0)
  }
  
  p <- p + ggplot2::geom_line()
  if (se) { ## add standard error if available
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se",
                                               ymax = ".S + .se",
                                               color = NULL,
                                               fill = color),
                           alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}