#' RIFplus
#'
#' Blah blah blah
#'
#' \tabular{ll}{ Package: \tab RIFplus\cr Type: \tab Package\cr Version:
#' \tab 0.1.0\cr Date: \tab 2016-12-19\cr License: \tab GPL (>=3)\cr LazyLoad:
#' \tab yes\cr }
#'
#' @name RIFplus-package
#' @aliases RIFplus-package
#' @docType package
#' @author Andrea Rau, Florence Jaffrezic, Toni Reverter
NULL


.onUnload <- function (libpath) {
  library.dynam.unload("RIFplus", libpath)
}
