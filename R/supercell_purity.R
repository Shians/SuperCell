#' Compute purity of super-cells
#'
#'
#' @param clusters vector of clustering assignment (reference assignment)
#' @param supercell_membership vector of assignment of single-cell data to super-cells (membership field of \code{\link{SCimplify}} function output)
#' @param ge gene expression use to compute the purity based on gene expression values
#'
#' @return a vector of super-cell purity, which is defined as a proportion of the most abundant cluster within super-cell.
#' With 1 meaning that super-cell consists of single cells from one cluster (reference assignment)
#'
#' @export
#'

supercell_purity <- function(ge = NULL, clusters=NULL, supercell_membership){
  if (!is.na(cluster)){
    cl.gr            <- table(clusters, supercell_membership)
    cluster.size     <- as.numeric(table(clusters))
    group.size       <- as.numeric(table(supercell_membership))

    Ng               <- length(group.size)
    group.max.cl     <- rep(0, Ng)


    cl.gr            <- sweep(cl.gr, 2, group.size, "/")

    res              <- apply(cl.gr, 2, max)
  } else {
    res <- rep(1, ncol(ge))
  }
  return(res)
}
