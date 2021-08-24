#' Compute purity of super-cells
#'
#'
#' @param clusters vector of clustering assignment (reference assignment)
#' @param supercell_membership vector of assignment of single-cell data to super-cells (membership field of \code{\link{SCimplify}} function output)
#'
#'
#' @return a vector of super-cell purity, which is defined as a proportion of the most abundant cluster within super-cell.
#' With 1 meaning that super-cell consists of single cells from one cluster (reference assignment)
#'
#' @export
#'

supercell_purity <- function(clusters, supercell_membership){
  cl.gr            <- table(clusters, supercell_membership)
  cluster.size     <- as.numeric(table(clusters))
  group.size       <- as.numeric(table(supercell_membership))

  Ng               <- length(group.size)
  group.max.cl     <- rep(0, Ng)


  cl.gr            <- sweep(cl.gr, 2, group.size, "/")

  res              <- apply(cl.gr, 2, max)
  return(res)
}


#' Purity based on transcriptionally
#'
#'  @param GE gene expression matrix
#'  @param membership membership vector
#'  @export

supercell_purity_GE <- function(GE, membership){
  # Compute cosine distance matrix
  ge <- sweep(GE, MARGIN = 2, STATS = sqrt(Matrix::colSums(GE**2)), FUN = "/")
  dist <- Matrix::crossprod(ge)

  # Make a membership mask matrix
  mask <- 1*(outer(membership, membership, FUN = "-") == 0)

  # Compute the mean purity
  purity <- Matrix::colSums(dist*mask)/Matrix::colSums(mask)
  return(purity)
}
