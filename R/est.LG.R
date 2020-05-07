#' Estimate graphons based on empirical degrees
#'
#' \code{est.LG} takes a 2-stage approach. First it adopts largest gap criterion on empirical degrees to
#' estimate blocks of a given network under Stochastic Blockmodel framework.
#' Then a consistent histogram estimator is utilized to estimate graphons based on
#' estimated blocks in a given network.
#'
#' @param A an \eqn{(n\times n)} binary adjacency matrix.
#' @param K the number of blocks provided by an user.
#'
#' @return a named list containing
#' \describe{
#' \item{H}{a  \eqn{(K\times K)} matrix of 3D histogram.}
#' \item{P}{an \eqn{(n\times n)} corresponding probability matrix.}
#' \item{B}{a length-\eqn{K} list where each element is a vector of nodes/indices
#' for each cluster.}
#' }
#'
#'
#' @examples
#' ## generate a graphon of type No.5 with 3 clusters
#' W = gmodel.preset(3,id=10)
#'
#' ## create a probability matrix for 20 nodes
#' graphW = gmodel.block(W,n=20)
#' P = graphW$P
#'
#' ## draw 23 observations from a given probability matrix
#' A = gmodel.P(P,rep=23,symmetric.out=TRUE)
#'
#' ## run LG algorithm with a rough guess for K=2,3,4
#' res2 = est.LG(A,K=2)
#' res3 = est.LG(A,K=3)
#' res4 = est.LG(A,K=4)
#'
#' ## compare true probability matrix and estimated ones
#' opar <- par(mfrow=c(2,2), pty="s")
#' image(P, main="original P matrix")
#' image(res2$P, main="LG with K=2")
#' image(res3$P, main="LG with K=3")
#' image(res4$P, main="LG with K=4")
#' par(opar)
#'
#' @references
#' \insertRef{Channarond2011}{graphon}
#'
#' \insertRef{chan2014}{graphon}
#'
#'
#' @seealso \code{\link{est.SBA}}
#' @export
est.LG <- function(A,K=2){
  ## Preprocessing : Directed Allowed
  if (is.vector(A)&&is.list(A)){
    if (!is.binAdjvec(A,sym=FALSE)){
      stop("* est.SBA : input matrix or vector A is invalid.")
    }
    vecA = A
  } else {
    if (!is.binAdj(A,sym=FALSE)){
      stop("* est.SBA : input matrix A is invalid.")
    }
    vecA = vector("list")
    vecA[[1]] = A
  }
  n = nrow(vecA[[1]])
  G = sum3(vecA,1:n,1:n,1:length(vecA))

  ## Preprocessing : K : number of blocks
  K = round(K)
  if ((K<1)||(K>n)){
    stop("* est.LG : the number of blocks K should be an integer in [1,number of nodes].")
  }

  ## Main Computation for Largest Gap Blocks
  Deg = rowSums(G-diag(diag(G)))           # degree
  DegNorm = Deg/(n-1)                      # normalized degree
  idx = order(DegNorm)                     # sorted index
  DegNorm = sort(DegNorm)                  # sort itself

  DegNormDiff = diff(DegNorm)              # difference
  Iod = order(DegNormDiff,decreasing=TRUE) # find (zk-1) largest gaps
  ii = c(0,sort(Iod[1:(K-1)]),n)

  B = list()
  for (k in 1:K){
    B[[k]] = idx[(ii[k]+1):ii[k+1]]
  }

  # histogram 3D
  result = histogram3D(vecA,B)
  result$B = B
  return(result)
}
