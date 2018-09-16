#' Generate graphs given a probability matrix
#'
#' Given an \eqn{(n\times n)} probability matrix \eqn{P}, \code{gmodel.P} generates
#' binary observation graphs corresponding to Bernoulli distribution
#' whose parameter matches to the element of \eqn{P}.
#'
#' @param P an \eqn{(n\times n)} probability matrix.
#' @param rep the number of observations to be generated.
#' @param noloop a logical value; TRUE for graphs without self-loops, FALSE otherwise.
#' @param symmetric.out a logical value; FALSE for generated graphs to be nonsymmetric, TRUE otherwise. Note that
#' TRUE is supported only if the input matrix P is symmetric.
#'
#' @examples
#' ## set inputs
#' modelP <- matrix(runif(16),nrow=4)
#'
#' ## generate 3 observations without self-loops.
#' out <- gmodel.P(modelP,rep=3,noloop=TRUE)
#'
#' ## Visualize generated graphs
#' par(mfrow=c(1,3))
#' image(out[[1]])
#' image(out[[2]])
#' image(out[[3]])
#'
#' @return depending on \code{rep} value, either
#' \describe{
#' \item{(rep=1)}{an \code{(n-by-n)} observation matrix, or}
#' \item{(rep>1)}{a length-\code{rep} list where each element
#' is an observation is an \code{(n-by-n)} realization from the model.}
#' }
#'
#' @export
gmodel.P <- function(P,rep=1,noloop=TRUE,symmetric.out=FALSE){
  ## Check P
  cond1 = ((all(P>=0))&&(all(P<=1)))
  cond2 = (nrow(P)==ncol(P))
  if (!(cond1&&cond2)){
    stop("* gmodel.P : P is not a valid probability matrix.")
  }

  ## Parameter
  n = nrow(P)

  ## Rep 1 case
  if (rep==1){
    tmpmat = matrix(runif(n^2),nrow=n)
    G = (tmpmat<P)*1
    if (noloop){
      diag(G) = 0
    }
  } else {
    G = list()
    for (i in 1:rep){
      tmpmat = matrix(runif(n^2),nrow=n)
      tmpG = 1*(tmpmat<P)
      if (noloop){
        diag(tmpG) = 0
      }
      G[[i]] = tmpG
    }
  }

  ## Symmetric
  if ((symmetric.out)&&(!isSymmetric(P))){
    stop("* gmodel.P : 'symmetric' option is only valid if where probability P matrix is symmetric.")
  }
  if (symmetric.out){
    if (isSymmetric(P)){
      if (rep==1){
        tmpG = (G+t(G))/2
        G = (tmpG>0)*1
      } else {
        for (i in 1:rep){
          tmptmpG = G[[i]]
          tmpG = tmptmpG+t(tmptmpG)
          G[[i]] = (tmpG>0)*1
        }
      }
    }
  }


  ## return output
  return(G)
}
