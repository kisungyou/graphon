#' Generate one of pre-specified graphons.
#'
#' \code{gmodel.preset} generates one of pre-specified graphons
#' of size \eqn{(n \times n)}. Users can select one of 10 different graphons by
#' their \code{id}, an integer from 1 to 10. The table of available graphons
#' follows that of the reference article given below.
#'
#' @param n the number of nodes for a graphon to be generated.
#' @param id an integer from 1 to 10, each corresponding to a specific graphon model.
#' @param sort a logical value; TRUE to sort in an decreasing order of degree, FALSE otherwise.
#'
#' @return an \eqn{(n\times n)} graphon matrix.
#'
#' @examples
#' \dontrun{
#' ## Generate 3 random graphons of nodal size 100.
#' n  = 100
#' r3 = (sample(1:10,3))
#' W1 = gmodel.preset(n,id=r3[1])
#' W2 = gmodel.preset(n,id=r3[2])
#' W3 = gmodel.preset(n,id=r3[3])
#'
#' ## Generate corresponding observations and plot them
#' A1 = gmodel.P(W1)
#' A2 = gmodel.P(W2)
#' A3 = gmodel.P(W3)
#'
#' \dontshow{
#'  for (i in 1:10){
#'    W = gmodel.preset(100,id=i)
#'  }
#' }
#' }
#'
#' @references
#' \insertRef{chan2014}{graphon}
#'
#' @export
gmodel.preset <- function(n, id=1, sort=TRUE){
  ## Parameters
  if (n<=1){
    stop("* gmodel.preset : the number of nodes should be >1.")
  }
  n = as.integer(n)
  if (!(id %in% 1:10)){
    stop("* gmodel.preset : graphon model id is an integer in [1,10].")
  }
  id = as.integer(id)

  ## Pass to CPP part
  vecgrid = seq(0,1,length.out=n)
  W = aux_preset(vecgrid,n,id)

  ## Decreasing degree
  if (sort){
    pos = order(colSums(W),decreasing=TRUE)
    Wtrue = W[pos,pos]
    return(Wtrue)
  } else {
    return(W)
  }
}
