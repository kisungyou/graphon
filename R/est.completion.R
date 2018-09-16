#' Estimate graphons based on matrix completion scheme
#'
#' \code{est.completion} adopts a matrix completion scheme,
#' which is common in missing data or matrix reconstruction studies.
#' When given a multiple of, or a single observation, we consider
#' non-existent edges as missing entries and apply the completion scheme.
#' See \code{\link[ROptSpace]{OptSpace}} for a more detailed introduction.
#'
#' @examples
#' ## generate a graphon of type No.5 with 3 clusters
#' W = gmodel.preset(3,id=5)
#'
#' ## create a probability matrix for 100 nodes
#' graphW = gmodel.block(W,n=100)
#' P = graphW$P
#'
#' ## draw 10 observations from a given probability matrix
#' A = gmodel.P(P,rep=10)
#'
#' ## apply the method
#' res_r3 = est.completion(A,rank=3)       # use rank-3 approximation
#' res_r9 = est.completion(A,rank=9)       # use rank-9 approximation
#' res_rN = est.completion(A,adjust=FALSE) # stop the code if guess works poorly
#'
#' ## visualize
#' par(mfrow=c(1,3))
#' image(res_r3); title("rank 3")
#' image(res_r9); title("rank 9")
#' image(res_rN); title("guessed rank")
#'
#' @references Keshavan, R.H., Montanari, A., and Oh, S. (2009) \emph{Matrix Completion
#' from a Few Entries}. Arxiv:0901.3150.
#'
#' @param A either \describe{
#' \item{Case 1.}{an \code{(n-by-n)} binary adjacency matrix, or}
#' \item{Case 2.}{a vector containing multiple of \code{(n-by-n)} binary adjacency matrices.}
#' }
#' @param rank an estimated rank condition for the matrix; \code{NA} for automatic guessing
#' of a rank, or a positive integer for a user-supplied rank assumption.
#' @param tolerance a tolerance level for singular value thresholding from OptSpace method.
#' @param maxiter the number of maximum iterations for OptSpace method.
#' @param progress a logical value; \code{FALSE} for not showing intermediate flags during
#' the process, \code{TRUE} otherwise.
#' @param adjust a logical value; \code{TRUE} to ignore a guessed rank and set it as 2 upon
#' numerical errors, \code{FALSE} to stop the code.
#'
#' @return an \code{(n-by-n)} corresponding probability matrix.
#'
#' @export
est.completion <- function(A,rank=NA,tolerance=1e-3,maxiter=20,progress=FALSE,adjust=TRUE){
  ## 1. Preprocessing : Directed Allowed
  if (is.vector(A)&&is.list(A)){
    if (!is.binAdjvec(A,sym=FALSE)){
      stop("* est.completion : input matrix or vector A is invalid.")
    }
    matA = array(0,c(nrow(A[[1]]),ncol(A[[2]])))
    for (i in 1:length(A)){
      matA = matA + A[[i]]
    }
    matA = matA/length(A)
  } else {
    if (!is.binAdj(A,sym=FALSE)){
      stop("* est.completion : input matrix A is invalid.")
    }
    matA = A
  }

  ## 2. Preprocessing : We consider zeroing out
  matA[(matA==0)] = NA
  OSinput = 2*(matA-0.5)

  ## 3. Run OptSpace
  optout = tryCatch({
    optout = OptSpace(OSinput,ropt=rank,tol=tolerance,niter=maxiter,showprogress=FALSE)
  }, warning = function(war){
    print("* Warning from OptSpace")
    print(paste("* Warning from OptSpace : ",war))
    if (adjust){
      print("* flag 'adjust' as TRUE : force the rank to be 2.")
      optout = OptSpace(OSinput,ropt=2,tol=tolerance,niter=maxiter,showprogress=progress)
    } else {
      print("* est.completion : stop for the message above.")
      stop()
    }
  }, error = function(war){
    print("* Error from OptSpace")
    print(paste("*",war))
    if (adjust){
      print("* flag 'adjust' as TRUE : force the rank to be 2.")
      optout = OptSpace(OSinput,ropt=2,tol=tolerance,niter=maxiter,showprogress=progress)
    } else {
      print("* est.completion : stop for the message above.")
      stop()
    }
  }
  )

  ## 4. get results
  X = optout$X
  S = optout$S
  Y = optout$Y

  output = ((X%*%S%*%t(Y))+1)/2
  output[(output>1)]=1
  output[(output<0)]=0
  return(output)
}
