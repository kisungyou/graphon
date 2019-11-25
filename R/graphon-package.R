#' graphon : A Collection of Graphon Estimation Methods
#'
#' The \pkg{graphon} provides a not-so-comprehensive list of methods for estimating graphon,
#' a symmetric measurable function, from a single or multiple of observed networks.
#' It also contains several auxiliary functions for generating sample networks using
#' various network models and graphons.
#'
#'
#' @section What is Graphon?:
#' Graphon - graph function - is a symmetric measurable function \deqn{W:[0,1]^2\rightarrow[0,1]} that arise
#' in studying exchangeable random graph models as well as sequence of dense graphs. In the language of
#' graph theory, it can be understood as a two-stage procedural network modeling that 1) each vertex/node in the graph
#' is assigned an independent random variable \eqn{u_j} from uniform distribution \eqn{U[0,1]}, and
#' 2) each edge \eqn{(i,j)} is randomly determined with probability \eqn{W(u_i,u_j)}. Due to such
#' procedural aspect, the term \emph{probability matrix} and \emph{graphon} will be interchangeably used
#' in further documentation.
#'
#'
#' @section Composition of the package:
#' The package mainly consists of two types of functions whose names start with \code{'est'}
#' and \code{'gmodel'} for estimation algorithms and graph models, respectively.
#'
#' The \code{'est'} family has 4 estimation methods at the current version,
#' \itemize{
#'   \item \code{\link{est.LG}} for empirical degree sorting in stochastic blockmodel.
#'   \item \code{\link{est.SBA}} for stochastic blockmodel approximation.
#'   \item \code{\link{est.USVT}} for universal singular value thresholding.
#'   \item \code{\link{est.nbdsmooth}} for neighborhood smoothing.
#'   \item \code{\link{est.completion}} for matrix completion from a partially revealed data.
#' }
#'
#' Also, the current release has following graph models implemented,
#' \itemize{
#'   \item \code{\link{gmodel.P}} generates a binary graph given an arbitrary probability matrix.
#'   \item \code{\link{gmodel.ER}} is an implementation of Erdos-Renyi random graph models.
#'   \item \code{\link{gmodel.block}} is used to generate networks with block structure.
#'   \item \code{\link{gmodel.preset}} has 10 exemplary graphon models for simulation.
#'  }
#'
#' @author Kisung You
#' @docType package
#' @name graphon-package
#' @import Rdpack
#' @importFrom utils packageVersion
#' @importFrom ROptSpace OptSpace
#' @importFrom graphics image par
#' @importFrom stats quantile rbinom runif
NULL



