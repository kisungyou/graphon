#' Generate binary random graphs based on stochastic blockmodel
#'
#' Given a \eqn{(K\times K)} stochastic blockmodel W, \code{gmodel.block}
#' generates an (n-by-n) binary random graphs. All K blocks have
#' same number of nodes, or almost identical if n is not a multiple
#' of K. Parameter \code{noloop} controls whether generated observations
#' have an edge from a node to itself, called a loop.
#'
#' @param W a \eqn{(K\times K)} blockmodel matrix.
#' @param n the number of nodes for each observation.
#' @param rep the number of observations to be generated.
#' @param noloop a logical value; TRUE for graphs without self-loops, FALSE otherwise.
#'
#' @seealso \code{\link{gmodel.P}}
#' @examples
#' ## set inputs
#' W = matrix(c(0.9,0.2,0.2,0.7),nr=2)
#' n = 200
#'
#' ## generate 2 observations without self-loops.
#' out <- gmodel.block(W,n,rep=2,noloop=TRUE)
#'
#' \donttest{
#' ## Visualize generated graphs
#' par(mfrow=c(1,2), pty="s")
#' image(out$G[[1]]); title("Observation 1")
#' image(out$G[[2]]); title("Observation 2")
#' }
#'
#' @return a named list containing
#' \describe{
#' \item{G}{depending on \code{rep} value, \describe{
#' \item{(rep=1)}{an \eqn{(n\times n)} observation, or}
#' \item{(rep>1)}{a length-\code{rep} list where each element
#' is an observation is an \eqn{(n\times n)} realization from the model.}
#' }
#' }
#' \item{P}{an \eqn{(n\times n)} probability matrix of generating each edge.}
#' }
#'
#' @export
gmodel.block <- function(W,n,rep=1,noloop=TRUE){
  ## Check W
  cond1 = ((all(W>=0))&&(all(W<=1)))
  cond2 = (nrow(W)==ncol(W))
  if (!(cond1&&cond2)){
    stop("* gmodel.block : W is not a valid stochastic block model.")
  }

  ## Parameter
  K = nrow(W) # number of blocks
  if (K>n){
    stop("* gmodel.block : the number of nodes should be >= the number of blocks.")
  }
  v = sort(runif(n))
  u = round((K-1)*v)+1

  ## Construct Probability Matrix
  P = array(0,c(n,n))
  for (i in 1:n){
    for (j in 1:n){
      ui = u[i]
      uj = u[j]
      P[i,j] = W[ui,uj]
    }
  }

  ## Construct a random graph with rep observations
  par_rep = rep
  par_loop = noloop
  G = gmodel.P(P,rep=par_rep,noloop=par_loop)

  ## output
  output = list()
  output$G = G
  output$P = P
  return(output)
}
