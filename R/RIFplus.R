#' Title
#'
#' @param x blah
#'
#' @return blah
#' @export
testfcn <- function(x) {
  out <- .Fortran("test", x = x, n = as.integer(length(x)),
                   PACKAGE = "RIFplus")
  return(list(out$x, out$n))
}


#' RIFplus function
#'
#' @param exprs Blah blah
#' @param TFnames Blah blah
#' @param conds Blah blah
#'
#' @return Value
#' @export
#' @importFrom methods is
#' @importFrom stats runif
#' @useDynLib RIFplus
RIFplus <- function(exprs, TFnames, conds) {

  ## Check: have data been normalized and log-transformed?
  ## Convert to data.frame, reorder columns, subset data
  if(is(exprs, "matrix")) exprs <- as.data.frame(exprs)
  conds <- factor(conds)
  o <- order(conds)
  conds <- conds[o]
  if(length(levels(conds))!=2)
    stop("Number of unique conditions must be = 2.")
  exprso <- exprs[,o]
  TFdata <- exprso[which(rownames(exprso) %in% TFnames),]
  targetdata <- exprso[which(!rownames(exprso) %in% TFnames),]
  ntf <- nrow(TFdata)
  nta <- nrow(targetdata)
  nconditions1 <- length(which(conds == levels(conds)[1]))
  nconditions2 <- length(which(conds == levels(conds)[2]))
  TFdata1 <- TFdata[,which(conds == levels(conds)[1])]
  TFdata2 <- TFdata[,which(conds == levels(conds)[2])]
  targetdata1 <- targetdata[,which(conds == levels(conds)[1])]
  targetdata2 <- targetdata[,which(conds == levels(conds)[2])]

  ## Convert data to vector (arranged by column)
  TFdata_vec1 <- as.numeric(unlist(t(TFdata1)))
  targetdata_vec1 <- as.numeric(unlist(t(targetdata1)))
  TFdata_vec2 <- as.numeric(unlist(t(TFdata2)))
  targetdata_vec2 <- as.numeric(unlist(t(targetdata2)))

  ## Initialize rif1scores and rif2scores
  rif1init <-  rif2init <- runif(ntf)
  out <- .Fortran("rif", ntf=as.integer(ntf),
                  nta=as.integer(nta),
                  nconditions1=as.integer(nconditions1),
                  nconditions2=as.integer(nconditions2),
                  TFdata1=TFdata_vec1,
                  targetdata1=targetdata_vec1,
                  TFdata2=TFdata_vec2,
                  targetdata2=targetdata_vec2,
                  rif1=rif1init,
                  rif2=rif2init,
                  cond1Xntf = as.integer(nconditions1 * ntf),
                  cond2Xntf = as.integer(nconditions2 * ntf),
                  cond1Xnta = as.integer(nconditions1 * nta),
                  cond2Xnta = as.integer(nconditions2 * nta),
                  PACKAGE = "RIFplus")
  rif1 <- round(out$rif1,5)
  rif2 <- round(out$rif2,5)
  rif1 <- (rif1 - mean(rif1)) /
    sqrt( (sum(rif1^2) - sum(rif1)*sum(rif1)/length(rif1)) / (length(rif1)-1) )
  rif2 <- (rif2 - mean(rif2)) /
    sqrt( (sum(rif2^2) - sum(rif2)*sum(rif2)/length(rif2)) / (length(rif2)-1) )
  return(list(
              rif1 = round(rif1,6), rif2 = round(rif2,6),
              rif1_top = data.frame(TFnames = TFnames[which(abs(rif1) >= 1.96 )],
                                    rif1 = rif1[which(abs(rif1) >= 1.96)],
                                    rif2 = rif2[which(abs(rif1) >= 1.96)]),
              rif2_top = data.frame(TFnames = TFnames[which(abs(rif2) >= 1.96 )],
                                    rif1 = rif1[which(abs(rif2) >= 1.96)],
                                    rif2 = rif2[which(abs(rif2) >= 1.96)]),
              rif1_nostd=out$rif1, rif2_nostd=out$rif2))
}




#' Example data
#'
#' Blah blah blah
#'
#' @name RIF_data
#' @docType data
#' @usage data(RIF_data)
#' @keywords datasets
#' @format Blah blah
#' @return Blah blah
NULL
