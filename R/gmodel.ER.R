#' Observations from Erdos-Renyi random graph model
#'
#' Erdos-Renyi random graph model is one of the most popular and
#' fundamental examples in modeling networks. Given n nodes,
#' \code{gmodel.ER} generates edges randomly from Bernoulli distribution
#' with a globally specified probability.
#'
#' In network science, 'ER' model is often interchangeably used in where
#' we have fixed number of edges to be placed at random. The original
#' use of edge-generating probability is from Gilbert (1959). Therefore,
#' we set this algorithm to be flexible in that user can create either a
#' fixed number of edges placed at random or set global edge-generating probability
#' and draw independent observations following Bernoulli distribution.
#'
#'
#' @param n the number of nodes to be generated
#' @param mode 'prob' (default) for edges to be drawn from Bernoulli
#' distribution independently, or 'num' for a graph to have a fixed
#' number of edges placed randomly
#' @param par a real number \eqn{\in [0,1]} for \code{mode='prob'}, or a
#' positive integer \eqn{\in [1, n*(n-1)/2]} for \code{mode='num'}
#' @param rep the number of observations to be generated.
#'
#' @return depending on \code{rep} value, either
#' \describe{
#' \item{(rep=1)}{an \eqn{(n\times n)} observation matrix, or}
#' \item{(rep>1)}{a length-\code{rep} list where each element
#' is an observation is an \eqn{(n\times n)} realization from the model.}
#' }
#'
#' @examples
#' ## generate 3 graphs with a global with probability 0.5
#' graph3 = gmodel.ER(100,par=0.5,rep=3)
#'
#' ## visualize
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' image(graph3[[1]], main="1st sample")
#' image(graph3[[2]], main="2nd sample")
#' image(graph3[[3]], main="3rd sample")
#' par(opar)
#'
#' @references
#' \insertRef{Erdos1959}{graphon}
#'
#' \insertRef{Gilbert1959}{graphon}
#'
#' @export
gmodel.ER <- function(n,mode='prob',par=0.5,rep=1){
  ## (1) Preprocessing
  ##  1-1. n
  if ((n<1)||is.na(n)||is.infinite(n)){
    stop('* gmodel.ER : number of nodes should be a positive integer.')
  }
  n = round(n)
  ##  1-2. mode detector
  if (!(mode %in% c('prob','num'))){
    stop('* gmodel.ER : mode should be either "prob" or "num."')
  }
  ##  1-3. mode-par
  maxpossible = round(n*(n-1)/2)
  if (mode=='prob'){
    if ((par<0)||(par>1)){
      stop('* gmodel.ER : for "prob" mode, par should be [0,1].')
    }
  } else{
    if ((par<1)){
      stop('* gmodel.ER : for "num" mode, par should be at least 1.')
    }
    if (par>=maxpossible){
      par = maxpossible
    }
    par = round(par)
  }
  ##  1-4. rep
  if ((rep<0)||(is.na(rep))||(is.infinite(rep))){
    stop('* gmodel.ER : number of replications should be a positive integer.')
  }
  rep = round(rep)

  ## (2) generation
  if (rep==1){
    if (mode=='prob'){
      A = slice.er.prob(maxpossible,par,n)
    } else {
      A = slice.er.num(maxpossible,par,n)
    }
  } else {
    A = vector("list",rep)
    for (i in 1:rep){
      if (mode=='prob'){
        tmpA = slice.er.prob(maxpossible,par,n)
      } else {
        tmpA = slice.er.num(maxpossible,par,n)
      }
      A[[i]] = tmpA
    }
  }

  ## (3) return output
  return(A)
}


# single slice generation for ER : probability method
slice.er.prob <- function(maxpossible,par,n){
  genvec = rbinom(maxpossible,1,par)
  A = matrix(0,n,n)
  A[upper.tri(A, diag=FALSE)] = genvec
  A = round(A+t(A))
  return(A)
}

# single slice generation for ER : #(edges) method
slice.er.num <- function(maxpossible,par,n){
  A = matrix(0,n,n)
  genvec = c(rep(1,par),rep(0,maxpossible-par))
  genvec[sample(1:maxpossible,maxpossible,replace=FALSE)] = genvec
  A[upper.tri(A, diag=FALSE)] = genvec
  A = round(A+t(A))
  return(A)
}
