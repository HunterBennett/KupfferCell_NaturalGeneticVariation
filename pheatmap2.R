library(pheatmap)

pheatmap2 <-
  function(D, cutree_rows = NA, ...) {
    #silent and annotation_row cannot be used
    if (!is.na(cutree_rows)) {
      X <- pheatmap(D, silent = T, ...)
      tr <- cutree(X$tree_row, cutree_rows)
      cluster <- paste("C", tr, sep = "")
      out <-
        pheatmap(
          D,
          annotation_row = data.frame(Clusters = cluster, row.names = rownames(D)),
          cutree_rows = cutree_rows,
          ...
        )
      cGene <- matrix("", ncol = 0, nrow = 0)
      for (i in 1:cutree_rows) {
        tmp <-
          matrix("",
                 nrow = max(nrow(cGene), sum(tr == i)),
                 ncol = ncol(cGene) + 1)
        if (nrow(cGene) > 0)
          tmp[1:nrow(cGene), 1:ncol(cGene)] <- cGene
        tmp[1:sum(tr == i), ncol(cGene) + 1] <- names(tr)[tr == i]
        cGene <- tmp
      }
      colnames(cGene) <- paste("C", 1:cutree_rows, sep = "")
      out <- c(out, list(clusterGene = cGene))
    } else{
      out <- pheatmap(D, ...)
    }
    return(out)
  }

pheatmap3 <-
  function(D, cutree_rows = NA, ...) {
    #silent and annotation_row cannot be used
    if (!is.na(cutree_rows)) {
      X <- pheatmap(D, silent = T, ...)
      tr <- cutree(X$tree_row, cutree_rows)
      cluster <- paste("C", tr, sep = "")
      out <-
        pheatmap(
          D,
          annotation_row = data.frame(Clusters = cluster, row.names = rownames(D)),
          cutree_rows = cutree_rows,
          ...
        )
    }
  }