#' Estimate edge probabilities by neighborhood smoothing
#'
#' \code{est.nbdsmooth} takes the expectation of the adjacency matrix
#' in that it directly aims at estimating network edge probabilities without
#' imposing structural assumptions as of usual graphon estimation requires,
#' such as piecewise lipschitz condition. Note that this method is
#' for symmetric adjacency matrix only, i.e., undirected networks.
#'
#'
#' @param A either \describe{
#' \item{Case 1.}{an \eqn{(n\times n)} binary adjacency matrix, or}
#' \item{Case 2.}{a vector containing multiple of \eqn{(n\times n)} binary adjacency matrices.}
#' }
#'
#' @return a named list containing
#' \describe{
#' \item{h}{a quantile threshold value.}
#' \item{P}{a matrix of estimated edge probabilities.}
#' }
#'
#' @examples
#' ## generate a graphon of type No.4 with 3 clusters
#' W = gmodel.preset(3,id=4)
#'
#' ## create a probability matrix for 100 nodes
#' graphW = gmodel.block(W,n=100)
#' P = graphW$P
#'
#' ## draw 5 observations from a given probability matrix
#' A = gmodel.P(P,rep=5,symmetric.out=TRUE)
#'
#' ## run nbdsmooth algorithm
#' res2 = est.nbdsmooth(A)
#'
#' ## compare true probability matrix and estimated ones
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' image(P, main="original P matrix")
#' image(res2$P, main="nbdsmooth estimated P")
#' par(opar)
#'
#' @references
#' \insertRef{Zhang2015}{graphon}
#'
#' @export
est.nbdsmooth <- function(A){
  # CASE 1. VECTORIAL CASE
  if (is.vector(A)&&is.list(A)){
    check1 = unlist(lapply(A,nrow))
    check2 = unlist(lapply(A,is.binAdj))

    cond1 = (length(unique(check1))==1)
    cond2 = (all(check2==TRUE))
    if (!(cond1&&cond2)){
      stop("* est.nbdsmooth : input vector A is not proper.")
    }
    # Main Computation for case 1
    nelem = nrow(A[[1]])
    P = matrix(0,nelem,nelem)

    for (i in 1:length(A)){
      tmpres = est.nbdsmoothsingle(A[[i]])
      P = P + tmpres$P
    }
    P = P/length(A)

    res = vector("list")
    res$h = tmpres$h
    res$P = P
  } else {
  # CASE 2. SINGLE CASE
    if (!is.binAdj(A)){
      stop("* est.nbdsmooth : input A is not a proper adjacency matrix.")
    }
    # Main Computation for case 2
    res = est.nbdsmoothsingle(A)
  }

  # Return Output
  return(res)
}


est.nbdsmoothsingle <- function(A){
  # 1. size
  N = nrow(A)
  h = sqrt(log(N)/N)

  # 2. compute dissimilarity measure
  D = aux_nbdsmooth(A, N)

  # 3. quantiled as logical
  kernel_mat = matrix(0,N,N)
  for (i in 1:N){
    kernel_mat[i,] = as.double(D[i,]<quantile(D[i,],h))
  }

  # 4. L1 normalization of each row
  kernel_mat = kernel_mat/(outer(rowSums(kernel_mat),rep(1,N))+1e-10)

  # 5. Compute P
  P = kernel_mat %*% A;
  P = (P+t(P))/2;

  ## (3) outputs
  res = vector("list")
  res$h = h
  res$P = P
  return(res)
}
