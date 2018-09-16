#' Cross validation for selecting optimal precision parameter in SBA method.
#'
#' The performance of Stochastic Blockmodel Approximation (SBA) method is
#' contingent on the number of blocks it finds during estimation process,
#' which is rougly determined by a precision parameter \code{delta}. \code{cv.SBA}
#' tests multiple of delta values to find the optimal one that minimizes
#' the cross validation risk. Note that the optimal delta is not bound to be a single value.
#'
#' @param A either \describe{
#' \item{Case 1.}{an \eqn{(n\times n)} binary adjacency matrix, or}
#' \item{Case 2.}{a vector containing multiple of \eqn{(n\times n)} binary adjacency matrices.}
#' }
#' @param vecdelta a vector containing target delta values to be tested.
#'
#' @return a named list containing \describe{
#' \item{optdelta}{optimal delta values that minimize the cross validation risk J.}
#' \item{J}{cross validation risk values.}
#' }
#'
#' @references
#' \insertRef{chan2014}{graphon}
#'
#' \insertRef{Airoldi2013}{graphon}
#'
#' @examples
#' \dontrun{
#' ## generate a graphon of type No.8 with 3 clusters
#' W = gmodel.preset(3,id=8)
#'
#' ## create a probability matrix for 100 nodes
#' graphW = gmodel.block(W,n=100)
#' P = graphW$P
#'
#' ## draw 15 observations from a given probability matrix
#' A = gmodel.P(P,rep=15)
#'
#' ## cross validate SBA algorithm over different deltas
#' rescv = cv.SBA(A,vecdelta=c(0.1,0.5,0.9))
#' print(rescv$optdelta)
#' }
#'
#'
#' @seealso \code{\link{est.SBA}}
#' @export
cv.SBA <- function(A,vecdelta=seq(0.1,1,by=0.1)){
  ## 1. Preprocessing : Directed Allowed
  if (is.vector(A)&&is.list(A)){
    if (!is.binAdjvec(A,sym=FALSE)){
      stop("* cv.SBA : input matrix or vector A is invalid.")
    }
    vecA = A
  } else {
    if (!is.binAdj(A,sym=FALSE)){
      stop("* cv.SBA : input matrix A is invalid.")
    }
    vecA = vector("list")
    vecA[[1]] = A
  }

  ## 2. Preprocessing : vecdelta
  ndelta = length(vecdelta)
  if (ndelta<=1){
    stop("* cv.SBA : vecdelta should be a vector of desired delta values for testing.")
  }

  ## 3.Main Iteration
  # 3-1. parameters and ready
  n = nrow(vecA[[1]])
  vecJ = array(0,c(ndelta,1))
  # 3-2. iterate over deltas
  for (i in 1:ndelta){
    # 3-2-1. estimate blocks
    cdelta = vecdelta[i]
    res = est.SBA(vecA,cdelta)
    B   = res$B

    # 3-2-2. compute p's
    K = length(B)
    if (K==1){
      vecp = 1
    } else {
      vecp = unlist(lapply(B,length))/n
    }

    # 3-2-3. compute J
    h = 1/K
    vecJ[i] = (2/(h*(n-1)))-(((n+1)/(h*(n-1)))*sum(vecp^2))
  }

  ## 4. results
  res = list()
  res$optdelta = vecdelta[(which(vecJ==min(vecJ)))]
  res$J = vecJ
  return(res)
}
