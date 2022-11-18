# XiaolinYan s2326461; Xiaoxian Zhang s2287323, Aysu Ismayilova s2295782
# https://github.com/AysuIsmayilova27/StatsProgramming4Group5/edit/main/proj4.R
# Everyone has contributed 1/3
newt <- function(theta,func,grad,hess=NULL,...,
                 tol=1e-8, fscale=1,maxit=100,max.half=20,eps=1e-6){
# write a newt function, using newton methods to get the minimization of function
# input:
# theta - a vector of initial values for the optimization parameters
# func - the objective function to minimize
# grad - the gradient function
# hess - the hessian matrix function
# tol - the convergence tolerance
# facale - a rough estimate of the magnitude of func near the optimum
# maxit - a vector of initial values for the optimization parameters
# max.half - the maximum number of times a step should be halved
# eps - the finite difference intervals to use to generate an approximate hessian
# output:
# f - the value of the objective function at the minimum
# theta - the value of the parameters at the minimum
# iter - the number of iterations taken to reach the minimum
# g - the gradient vector at the minimum
# Hi - the inverse of the Hessian matrix at the minimum 

  
  f0 <- func(theta,...) # the initial value of our objective function
  gradient <- grad(theta,...) # the initial gradient of intial objective function
  
  
  if(is.null(hess)){ # if hessian matrix is not provided,we need to generate an approximation
    hess <- function(theta) {
      eps <- 1e-8 # finite difference intervals
      gradient <- grad(theta)
      fd.f <- matrix(0, nrow = length(theta), ncol = length(theta)) # initialize hessian matrix
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
    stop("The objective or its derivatives are not finite at initial theta.")
  }
  
  # optimization begins
  iter = 0 # initialize iteration
  while (iter < maxit) {
    if (max(abs(gradient)) < (abs(f0)+fscale)*tol){# if converged
      cat("Converged",fill=TRUE)
      # check if hessian matrix is positive definte at convergence,we cannot get cholesky decomposition if hessian is not positive definite
      if (inherits(try(chol(H), silent = TRUE),"try-error")){ # if try error, hessian matrix is not positive definite, give warning
        warning("The Hessian is not positive definite at convergence")
      }
     
      Hi <- ginv(H) # inverse hessian matrix
      # return a list contain the minimum f,the value of theta, iteration, the gradient and inverse of the Hessian matrix at the minimum
      return(list(f=f0[1], theta=theta, iter=iter, g=gradient, Hi=Hi))
    } else { # if not converged
      catch.iter <- 0
      while(inherits(try(cholesky<- chol(H), silent = TRUE), "try-error") == TRUE) { 
       #if try error, hessian matrix is not positive definite, we should make it positive definite
        H <- H + diag(abs(max(H))*10^catch.iter*tol, nrow=nrow(H), ncol=ncol(H))
        catch.iter <- catch.iter + 1
      }
      # newton step, we need to calculate forward step to optimize our function
      Delta <- backsolve(cholesky, forwardsolve(t(cholesky), -gradient))
     
      # judge if the step fails to reduce objective function after trying max.half step halving
      for (i in 1:max.half){
        f1 <- func(theta+Delta,...)
        g1 <- grad(theta + Delta,...)
        h1 <- hess(theta + Delta)
        
        if((f1[1] > f0[1]) |is.infinite(f1) |any(is.infinite(g1)) |any(is.infinite(h1))){
          Delta <- Delta / 2 # step halve
        }
      }
      
      theta <- theta + Delta # Updating theta
      f0 <- func(theta,...) #  Updating objective function 
      gradient <- grad(theta,...) #  Updating gradient
      H <- hess(theta)#  Updating hessian matrix
      iter <- iter + 1 # updating iter
      # give warning If the step fails to reduce the objective despite trying max.half step halvings
      if(iter==max.half){
        warning(" fails to reduce the objective despite trying max.half step halvings")
      }
    }
  }
  # iteration limit
 if(iter==maxit){
    if (max(abs(gradient)) < (abs(f0)+fscale)*tol){# If gradient is very close to the limit
      cat("Converged",fill=TRUE)# It converged
      return(list(f= f0, theta =theta, iter =iter, g =gradient, Hi= Hi))# And we return the list 
    } else {
      warning("Newton optimizer failed to converge after maxit = ", as.character(maxit), " iterations")
      # reach the maximum iterations, we give a warning
    }
  }
  
  
}
