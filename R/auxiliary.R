#' @keywords internal
#' @noRd
is.Adj <- function(A){
  # 1. size
  cond1 = ((is.matrix(A))&&(nrow(A)==ncol(A)))
  # 2. symmetric : bounded
  cond2 = (sum((A-t(A))^2)<1e-10)
  # 3. no negative values
  cond3 = (all((A>=0)))
  # 4. diagonals are zeros
  cond4 = (all(diag(A)==0))

  if (cond1&&cond2&&cond3&&cond4){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

## Binary {0,1} adjacency
#' @keywords internal
#' @noRd
is.binAdj <- function(A,sym=TRUE){
  # 1. size
  cond1 = ((is.matrix(A))&&(nrow(A)==ncol(A)))
  # 2. symmetric : bounded
  if (sym){
    cond2 = (sum((A-t(A))^2)<1e-10)
  } else {
    cond2 = TRUE
  }
  # 3. no negative values
  cond3 = (all((A>=0)))
  # 4. diagonals are zeros
  cond4 = (all(diag(A)==0))
  # 5. all binaries
  cond5 = all((A==1)||(A==0))

  if (cond1&&cond2&&cond3&&cond4&&cond5){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' @keywords internal
#' @noRd
is.binAdjvec <- function(vecA,sym=TRUE){
  # 1. size
  cond1 = (length(unique(unlist(lapply(vecA,nrow)))==1))
  # 2. other factors
  if (sym){
    symvec = TRUE
  } else {
    symvec = FALSE
  }

  cvec = array(0,c(length(vecA),1))
  for (i in 1:length(vecA)){
    if (is.binAdj(vecA[[i]],sym=symvec)){
      cvec[i] = 1
    }
  }
  if (sum(cvec)==length(cvec)){
    cond2 = TRUE
  } else {
    cond2 = FALSE
  }

  if (cond1&&cond2){
    return(TRUE)
  } else{
    return(FALSE)
  }
}

#' @keywords internal
#' @noRd
sum3 <- function(vecA,x,y,z){
  T = length(vecA)
  if (is.logical(x)){
    nx = which(x)
  } else {
    nx = length(x)
  }
  if (is.logical(y)){
    ny = which(y)
  } else {
    ny = length(y)
  }
  res = array(0,c(nx,ny))
  for (i in 1:length(z)){
    tgtnum = z[i]
    res = res+vecA[[tgtnum]][x,y]
  }
  return(res)
}

# Inputs
#   G : vector list of (n x n) graph
#   B : vector list of clusters
# Output - a list containing
#   H : 3D histogram
#   P : corresponding probability matrix of (n x n)
#' @keywords internal
#' @noRd
histogram3D <- function(G,B){
  # 1. get information about data
  n  = nrow(G[[1]])
  nT = length(G)
  nB = length(B)

  # 2. initialization
  P = array(0,c(n,n))
  H = array(0,c(nB,nB))

  # 3. Loop through all the clusters
  for (ki in 1:nB){
    for (kj in 1:nB){
      # 3-1. obtain the indices in cluster ki and kj
      I = B[[ki]]
      J = B[[kj]]

      # 3-2. compute 3D histogram
      H[ki,kj] = sum(sum3(G,I,J,1:nT))/(T*length(I)*length(J))

      # 3-3. loop through the indices in cluster I and J
      # to compute the corresponding probability matrix
      for (i1 in 1:length(I)){
        for (j1 in 1:length(J)){
          vi = I[i1]
          vj = J[j1]
          P[vi,vj] = H[ki,kj]
        }
      }
    }
  }

  # 4. return output in a list format
  res = list()
  res$H = H
  res$P = P
  return(res)
}

#' @keywords internal
#' @noRd
aux_nbdsmooth <- function(A,N){
  # definitions
  D = array(0,c(N,N))
  A_sq = (A%*%A)/N

  # compute D : dissimilarity
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      # val = max(abs(A_sq[i,]-A_sq[j,])) # this is my original work
      tgtvec = abs(A_sq[i,]-A_sq[j,])
      tgtvec[i] = 0
      tgtvec[j] = 0
      val = max(tgtvec) # tgtvec2 is Li Chen's
      D[i,j] = val
      D[j,i] = val
    }
  }

  # return result
  return(D)
}

# gmodel.preset : fast generation -----------------------------------------
#' @keywords internal
#' @noRd
aux_preset <- function(vecgrid,n,id){
  # preparation
  W = array(0,c(n,n))

  # iteration
  for (i in 1:n){
    u = vecgrid[i]
    for (j in 1:n){
      v = vecgrid[j]
      if (id==1){
        W[i,j] = u*v;
      } else if (id==2){
        W[i,j] = exp(-((u**0.7)+(v**0.7)));
      } else if (id==3){
        W[i,j] = ((u**2)+(v**2)+sqrt(u)+sqrt(v))/4;
      } else if (id==4){
        W[i,j] = (u+v)/2;
      } else if (id==5){
        W[i,j] = 1/(1+exp(-10*((u**2)+(v**2))));
      } else if (id==6){
        W[i,j] = abs(u-v)
      } else {
        maxuv = max(u,v)
        minuv = min(u,v)
        if (id==7){
          W[i,j] = 1/(1+exp(-((maxuv**2)+(minuv**4))));
        } else if (id==8){
          W[i,j] = exp(-(maxuv**0.75));
        } else if (id==9){
          W[i,j] = exp(-0.5*(minuv+sqrt(u)+sqrt(v)));
        } else if (id==10){
          W[i,j] = log(1+0.5*maxuv);
        }
      }
    }
  }

  # return W
  return(W)
}
