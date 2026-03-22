#' @keywords internal
"_PACKAGE"

#' @useDynLib hhdynamics, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom stats model.matrix quantile reshape sd density median acf pweibull qnorm
#' @importFrom graphics par plot abline polygon lines segments points barplot arrows axis legend
#' @importFrom grDevices adjustcolor grey.colors pdf dev.off
NULL
