# All code credit to Corey Chivers (github: cjbayesian) and the MHadaptive Library 
# https://github.com/cjbayesian/MHadaptive/blob/master/R/positiveDefinite.R

# functions I need for adaptive mcmc

.make.positive.definite <-
  function(m, tol)
  {
    if (!is.matrix(m)) {
      m = as.matrix(m)
    }
    
    d = dim(m)[1]
    if ( dim(m)[2] != d ) {
      stop("Input matrix is not square!")
    }
    
    es = eigen(m)
    esv = es$values
    
    if (missing(tol)) {
      tol = d*max(abs(esv))*.Machine$double.eps
    }
    delta =  2*tol
    
    tau = pmax(0, delta - esv)
    dm = es$vectors %*% diag(tau, d) %*% t(es$vectors)
    
    return( m + dm )
  }


makePositiveDefinite <-
  function(x)
  {
    # Arguments:
    #   x - a symmetric matrix or any other rectangular object
    #       describing a covariance matrix which can be converted into
    #       a matrix by the function 'as.matrix'.
    
    ans = .make.positive.definite(m = x)
    
    ans
  }