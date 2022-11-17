
newt <- function(theta,func,grad,hess=NULL,...,
                 tol=1e-8, fscale=1,maxit=100,max.half=20,eps=1e-6){
  
  f0 <- func(theta,...) # the initial value of our objective function
  gradient <- grad(theta,...) # the initial gradient of intial objective function
  
  
  if(is.null(hess)){ # if hessian matrix is not provided,we need to generate an approximation
    hess <- function(theta) {
      eps <- 1e-8 # finite difference intervals
      gradient <- grad(theta)
      fd.f <- matrix(0L, nrow = length(theta), ncol = length(theta)) # initialize hessian matrix
      for (i in 1:length(theta)){ 
        th1 <- theta; th1[i] <- th1[i]+eps
        f.hi <- grad(th1)
        fd.f[,i] <- (f.hi - gradient)/eps # generate an approximate hessian matrix by finite differencing of the gradient vector 
      }
      
      fd.f <- 0.5 * (t(fd.f) + fd.f) # get symmetric approximate hessian matrix
      return(fd.f)
    }
    H <- hess(theta)

  }else{ # if hessian matrix is provided, get the initial hessian matrix of intial objective function
  H <- hess(theta,...)
  H <- 0.5 * (t(H) + H) # make hessian matrix symmetric
  }
  
  #  If the objective or its derivatives are not finite at the initial theta, stop
  if (is.finite(f0)==FALSE | any(is.finite(gradient)==FALSE)| any(is.finite(H)==FALSE)) {
    stop("The objective or its derivatives are not finite at
         your initial estimate of theta.")
  }
  
  # optimization begins
  iter = 0 # initialize iteration
  while (iter < maxit) {
    if (max(abs(gradient)) < (abs(f0)+fscale)*tol){# if converged
      cat("Converged",'\n')
      # check if hessian matrix is positive definte at convergence,we cannot get cholesky decomposition if hessian is not positive definite
      if (inherits(try(chol(H), silent = TRUE),"try-error")){ # if try error, hessian matrix is not positive definite, give warning
        warning("The Hessian is not positive definite at convergence")
      }
     
      Hi <- ginv(H) # inverse hessian matrix
      # return a list contain the minimum f,the value of theta, iteration, the gradient and inverse of the Hessian matrix at the minimum
      return(list(f.min=f0[1], theta=theta, iter=iter, g=gradient, Hi=Hi))
    } else { # if not converged
      catch.iter <- 0
      while(inherits(try(H.chol <- chol(H), silent = TRUE), "try-error") == TRUE) { 
       #if try error, hessian matrix is not positive definite, we should make it positive definite
        H <- H + diag(abs(max(H))*10^catch.iter*tol, nrow=nrow(H), ncol=ncol(H))
        catch.iter <- catch.iter + 1
      }
      # newton step, we need to calculate forward step to optimize our function
      Delta <- backsolve(H.chol, forwardsolve(t(H.chol), -gradient))
     
      # judge if the step fails to reduce objective function after trying max.half step halving
      half.iter = 0 # initialize interation
     
      while ((func(theta + Delta)[1] > f0[1]) |
             is.infinite(func(theta + Delta)) |
             any(is.infinite(grad(theta + Delta))) |
             any(is.infinite(hess(theta + Delta)))) {
        if (half.iter < max.half) {
          Delta <- Delta / 2
          half.iter <- half.iter + 1
        } else {
          stop(paste("The update step failed to reduce the objective
                  after ", as.character(max.half), " halvings"))
        }
      }

      theta <- theta + Delta
      
     
      f0 <- func(theta)
      gradient <- grad(theta)
      
      
      if (any(is.null(hess(theta,...)))) {
       
        H <- fd.f
      } else {
        H <- hess(theta,...)
      }
      

      iter <- iter + 1
    }
  }
 
  if (max(abs(gradient)) < (abs(f0)+fscale)*tol){
    cat("Converged")
    return(list(f0, theta, iter, gradient, Hi))
  } else {
    warning(paste("Newton optimizer failed to converage after
                  maxit = ", as.character(maxit), " iterations"))
  }
  
  
}
