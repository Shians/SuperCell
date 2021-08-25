#' Compute the sparsity of gene counts
#'
#' @param GE
#' @param genes the set of gene to inspect
#' @export

supercell_NonZero <- function(GE, membership, genes){
  coverage_sc <- Matrix::rowSums(1*(GE[genes, ] != 0))
  GE_SC <- supercell_GE(GE, membership)
  coverage_SC <- Matrix::rowSums(1*(GE_SC[genes, ] != 0))/dim(GE_SC)[2]
  return(coverage_SC)
}
