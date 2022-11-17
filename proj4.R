





library('MASS')

newt <- function(theta,func,grad,hess=NULL,...,
                 tol=1e-8, fscale=1,maxit=100,max.half=20,eps=1e-6){
  
  f0 <- func(theta,...)
  gradient <- grad(theta,...)
  if(is.null(hess)){
    hess <- function(theta) {
      eps <- 1e-8
      gradient <- grad(theta)
     
      fd.f <- matrix(0L, nrow = length(theta), ncol = length(theta))
      for (i in 1:length(theta)){
        th1 <- theta; th1[i] <- th1[i]+eps
        f.hi <- grad(th1)
        fd.f[,i] <- (f.hi - gradient)/eps
      }
      
      fd.f <- 0.5 * (t(fd.f) + fd.f)
      return(fd.f)
    }
    H <- hess(theta)

  }else{
  H <- hess(theta,...)
  H <- 0.5 * (t(H) + H)
  }

  if (is.finite(f0)==FALSE | any(is.finite(gradient)==FALSE)| any(is.finite(H)==FALSE)) {
    stop("The objective or its derivatives are not finite at
         your initial estimate of theta.")
  }

  iter = 0
  while (iter < maxit) {
    if (max(abs(gradient)) < (abs(f0)+fscale)*tol){# Checking convergence
      cat("Converged",'\n')
      if (inherits(try(chol(H), silent = TRUE),"try-error")){ 
       
        warning("The Hessian is not positive definite at convergence")
      }
     
      Hi <- ginv(H)
      return(list(f.min=f0[1], theta=theta, iter=iter, g=gradient, Hi=Hi))
    } else {
      catch.iter <- 0
     
      while(inherits(try(H.chol <- chol(H), silent = TRUE), "try-error") == TRUE) {
        
        H <- H + diag(abs(max(H))*10^catch.iter*tol, nrow=nrow(H), ncol=ncol(H))
        catch.iter <- catch.iter + 1
      }
     
      Delta <- backsolve(H.chol, forwardsolve(t(H.chol), -gradient))
      
     
      half.iter = 0
      
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
th0 <- c(0.5,2)

rb <- function(th,k=2) {# func = objective function
  # D = k*(x_2 - x_1^2)^2 + (1-x_1)^2
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2 
}

gb <- function(th,k=2) {# grad = gardient
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}


hb <- function(th,k=2) {# hess = Hessian matrix
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}

newt(th0,rb,gb,hb) ## full version
