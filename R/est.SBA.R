#' Estimate graphons based on Stochastic Blockmodel approximation
#'
#' \code{est.SBA} takes a 2-stage approach for estimating graphons
#' based on exchangeable random graph models. First, it finds a
#' Stochastic Blockmodel Approximation (SBA) of the graphon. Then,
#' it uses clustering information to estimate graphon using a consistent
#'  histogram estimator.
#'
#' @param A either \describe{
#' \item{Case 1.}{an \eqn{(n\times n)} binary adjacency matrix, or}
#' \item{Case 2.}{a vector containing multiple of \eqn{(n\times n)} binary adjacency matrices.}
#' }
#' @param delta a precision parameter larger than 0.
#'
#' @return a named list containing
#' \describe{
#' \item{H}{a \eqn{(K\times K)} matrix fo 3D histogram.}
#' \item{P}{an \eqn{(n\times n)} corresponding probability matrix.}
#' \item{B}{a length-\eqn{K} list where each element is a vector of nodes/indices
#' for each cluster.}
#' }
#'
#' @examples
#' \donttest{
#' ## generate a graphon of type No.6 with 3 clusters
#' W = gmodel.preset(3,id=6)
#'
#' ## create a probability matrix for 100 nodes
#' graphW = gmodel.block(W,n=100)
#' P = graphW$P
#'
#' ## draw 17 observations from a given probability matrix
#' A = gmodel.P(P,rep=17)
#'
#' ## run SBA algorithm with different deltas (0.2,0.5,0.8)
#' res2 = est.SBA(A,delta=0.2)
#' res3 = est.SBA(A,delta=0.5)
#' res4 = est.SBA(A,delta=0.8)
#'
#' ## compare true probability matrix and estimated ones
#' par(mfrow=c(1,4))
#' image(P); title("original P")
#' image(res2$P); title("SBA with delta=0.2")
#' image(res3$P); title("SBA with delta=0.5")
#' image(res4$P); title("SBA with delta=0.8")
#' }
#'
#' @references
#' \insertRef{Airoldi2013}{graphon}
#'
#' \insertRef{chan2014}{graphon}
#'
#'
#' @seealso \code{\link{est.LG}}
#' @export
est.SBA <- function(A,delta=0.5){
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
  ## Preprocessing : delta
  if (delta<=0){
    stop("* est.SBA : delta should be a positive number.")
  }

  ## Main Estimation
  # 1. parameters
  nT = length(vecA)    # number of random graphs
  n  = nrow(vecA[[1]]) # number of nodes
  B  = vector("list")  # empty block labels

  PivotIdx = sample(1:n,1) # initial pivot index
  B[[1]]   = PivotIdx      # Block 1 should contain Pivot 1

  NotAssignedVector = rep(TRUE,n)
  NotAssignedVector[PivotIdx] = FALSE
  NotAssignedIdx = which(NotAssignedVector)

  # 2. main iteration
  #   loops until
  #   - all indices have been assigned
  #   - all nodes have been scanned
  while (any(NotAssignedVector==TRUE)){
    # 2-1. randomly choose a vectex from Omega
    if (length(NotAssignedIdx)>1){
      i = sample(NotAssignedIdx,1)
    } else {
      i = NotAssignedIdx
    }
    NotAssignedVector[i] = FALSE
    NotAssignedIdx = which(NotAssignedVector)

    # 2-2. Compute the distance estimate d
    dhat = array(0,c(length(PivotIdx),1))
    for (j in 1:length(PivotIdx)){
      # define the j-th pivot
      bj = PivotIdx[j]

      # define the set S (nbd for computing dhat)
      SVector = (!logical(length=n))
      SVector[c(i,bj)] = FALSE
      SIdx = which(SVector)

      # compute dhat using formula (5) in SBA paper
      Tf = floor((nT+1)/2)
      if (Tf<nT){
        Term1 = sum(((1/Tf)*sum3(vecA, i,SIdx,1:Tf))*((1/(nT-Tf))*sum3(vecA, i,SIdx,(Tf+1):nT)))
        Term2 = sum(((1/Tf)*sum3(vecA,bj,SIdx,1:Tf))*((1/(nT-Tf))*sum3(vecA,bj,SIdx,(Tf+1):nT)))
        Term3 = sum(((1/Tf)*sum3(vecA, i,SIdx,1:Tf))*((1/(nT-Tf))*sum3(vecA,bj,SIdx,(Tf+1):nT)))
        Term4 = sum(((1/Tf)*sum3(vecA,bj,SIdx,1:Tf))*((1/(nT-Tf))*sum3(vecA,i, SIdx,(Tf+1):nT)))

        Term5 = sum(((1/Tf)*sum3(vecA,SIdx, i,1:Tf))*((1/(nT-Tf))*sum3(vecA,SIdx,i, (Tf+1):nT)))
        Term6 = sum(((1/Tf)*sum3(vecA,SIdx,bj,1:Tf))*((1/(nT-Tf))*sum3(vecA,SIdx,bj,(Tf+1):nT)))
        Term7 = sum(((1/Tf)*sum3(vecA,SIdx, i,1:Tf))*((1/(nT-Tf))*sum3(vecA,SIdx,bj,(Tf+1):nT)))
        Term8 = sum(((1/Tf)*sum3(vecA,SIdx,bj,1:Tf))*((1/(nT-Tf))*sum3(vecA,SIdx,i, (Tf+1):nT)))
      } else {
        Term1 = sum(((1/Tf)*sum3(vecA, i,SIdx,1:Tf)))
        Term2 = sum(((1/Tf)*sum3(vecA,bj,SIdx,1:Tf)))
        Term3 = sum(((1/Tf)*sum3(vecA, i,SIdx,1:Tf)))
        Term4 = sum(((1/Tf)*sum3(vecA,bj,SIdx,1:Tf)))

        Term5 = sum(((1/Tf)*sum3(vecA,SIdx, i,1:Tf)))
        Term6 = sum(((1/Tf)*sum3(vecA,SIdx,bj,1:Tf)))
        Term7 = sum(((1/Tf)*sum3(vecA,SIdx, i,1:Tf)))
        Term8 = sum(((1/Tf)*sum3(vecA,SIdx,bj,1:Tf)))
      }

      dhatTmp = 0.5*(abs(Term1+Term2-Term3-Term4) + abs(Term5+Term6-Term7-Term8));
      dhat[j] = sqrt(abs(dhatTmp/length(SIdx)));
    }

    # 2-3. Assign Clusters
    Val = min(dhat)
    Idx = which(dhat==Val)
    if (Val<delta){
      # If min distance < Delta, assign to one of the existing blocks
      B[[Idx]] = c(B[[Idx]], i)
    } else {
      # If min distance > Delta, make a new block; Put i as pivot
      B[[length(PivotIdx)+1]] = i
      PivotIdx = c(PivotIdx, i)
    }
  }

  # 3. histogram 3D
  result = histogram3D(vecA,B)
  result$B = B
  return(result)
}
