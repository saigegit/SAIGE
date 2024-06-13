# l2.var: the variance of the prior, so the penalty is 1/(2*l2.var)
# l2.var.pos: the position of the covariates that applied to l2 penalty
fast.logistf.fit.WithL2 <- function (x, y, weight = NULL, offset = NULL, firth = TRUE, col.fit = NULL,
                              init = NULL, l2.var=NULL, l2.var.pos=NULL, control) {
  n <- nrow(x)
  k <- ncol(x)
  if (is.null(init))
    init = rep(0, k)
  if (is.null(col.fit))
    col.fit = 1:k
  if (is.null(offset))
    offset = rep(0, n)
  if (is.null(weight))
    weight = rep(1, n)
  if (col.fit[1] == 0)
    maxit <- 0
  if (missing(control))
    control <- fast.logistf.control()
  maxit <- control$maxit
  maxstep <- control$maxstep
  maxhs <- control$maxhs
  lconv <- control$lconv
  gconv <- control$gconv
  xconv <- control$xconv
  beta <- init
  iter <- 0
  pi <- as.vector(1/(1 + exp(-x %*% beta - offset)))
  evals <- 1
  repeat {
    beta.old <- beta
    XW2 <- t(x * (weight * pi * (1-pi))^0.5)
    myQR <- qr(t(XW2))
    Q <- qr.Q(myQR)
    h <- (Q*Q) %*% rep(1, ncol(Q))
    if (firth)
      U.star <- crossprod(x, weight * (y - pi) + h * (0.5 - pi))
    else U.star <- crossprod(x, weight * (y - pi))

    # Added by SLEE 2023/07/22
    if(!is.null(l2.var)){
      U.star[l2.var.pos] = U.star[l2.var.pos] - 1/l2.var *beta[l2.var.pos]
    }
    ##

    XX.covs <- matrix(0, k, k)

    if (col.fit[1] != 0) {
      XX.XW2 <- t(x[, col.fit, drop=FALSE] * (weight * pi * (1-pi))^0.5)
      XX.Fisher <- crossprod(t(XX.XW2))

      # Added by SLEE 2023/07/22
      if(!is.null(l2.var)){
        XX.Fisher[l2.var.pos, l2.var.pos] = XX.Fisher[l2.var.pos, l2.var.pos] +  1/l2.var
      }
      ##

      XX.covs[col.fit, col.fit] <- fast.invFisher(XX.Fisher)  ###### HERE IS THE PROBLEM!!!
    }
    if(all(is.na(XX.covs)) == T) {
      break
    }
    delta <- as.vector(XX.covs %*% U.star)
    delta[is.na(delta)] <- 0
    mx <- max(abs(delta))/maxstep
    if (mx > 1)
      delta <- delta/mx
    evals <- evals + 1
    if (maxit > 0) {
      iter <- iter + 1
      beta <- beta + delta
      pi <- as.vector(1/(1 + exp(-x %*% beta - offset)))

    }
    if (iter == maxit | ((max(abs(delta)) <= xconv) & (all(abs(U.star[col.fit]) <
                                                           gconv))))
      break
  }
  # Error catching (if chol(x) not positive definite)
  if(all(is.na(XX.covs))==T) {
    var <- XX.covs
    list(beta = NA, var = var, pi = NA, hat.diag = NA,
         iter = NA, evals = NA, conv = c(NA,
                                         NA, NA))
  } else {
    var <- XX.covs
    list(beta = beta, var = var, pi = pi, hat.diag = h,
         iter = iter, evals = evals, conv = c(max(abs(U.star)),
                                              max(abs(delta))))
  }
}


# By SLEE
# calculate h under the null, so intercept is
calculate_firth_null <- function (y, weight = NULL, offset = NULL) {

  n <- length(y)
  if (is.null(offset))
    offset = rep(0, n)
  if (is.null(weight))
    weight = rep(1, n)

  out.null.firth = fast.logistf.fit.WithL2(cbind(rep(1,n)), y, weight=weight, offset=offset)
  pi <- out.null.firth$pi
  out.null.firth$pi_var_sum = sum(pi* (1-pi))
 return(out.null.firth)

}


# By SLEE
# Do not update fisher information
# l2.var: the variance of the prior, so the penalty is 1/(2*l2.var)
# l2.var.pos: the position of the covariates that applied to l2 penalty
fast.logistf.fit.WithL2.fast <- function (x, y, out.null.firth=NULL, weight = NULL, offset = NULL, firth = TRUE, col.fit = NULL,
                                     init = NULL, l2.var=NULL, l2.var.pos=NULL, control) {
  n <- nrow(x)
  k <- ncol(x)
  if (is.null(init))
    init = rep(0, k)
  if (is.null(col.fit))
    col.fit = 1:k
  if (is.null(offset))
    offset = rep(0, n)
  if (is.null(weight))
    weight = rep(1, n)
  if (col.fit[1] == 0)
    maxit <- 0
  if (missing(control))
    control <- fast.logistf.control()
  maxit <- control$maxit
  maxstep <- control$maxstep
  maxhs <- control$maxhs
  lconv <- control$lconv
  gconv <- control$gconv
  xconv <- control$xconv

  re.out<-list(beta = NA, var = var, pi = NA, hat.diag = NA, iter = NA, evals = NA, conv = c(NA, NA, NA))

  # By SLEE, compute once...
  if(is.null(out.null.firth)){
  	out.null.firth = fast.logistf.fit.WithL2(cbind(rep(1,n)), y, weight=weight, offset=offset)
  }

  pi <- out.null.firth$pi
  beta <- c(out.null.firth$beta,0)
  iter <- 0
  evals <- 1


  XW2 <- t(x * (weight * pi * (1-pi))^0.5)
  myQR <- qr(t(XW2))
  Q <- qr.Q(myQR)
  h <- (Q*Q) %*% rep(1, ncol(Q))


  repeat {

    beta.old <- beta
    if (firth)
      U.star <- crossprod(x, weight * (y - pi) + h * (0.5 - pi))
    else U.star <- crossprod(x, weight * (y - pi))

    # Added by SLEE 2023/07/22
    if(!is.null(l2.var)){
      U.star[l2.var.pos] = U.star[l2.var.pos] - 1/l2.var *beta[l2.var.pos]
    }
    ##
    XX.covs <- matrix(0, k, k)
    if (col.fit[1] != 0) {
      XX.XW2 <- t(x[, col.fit, drop=FALSE] * (weight * pi * (1-pi))^0.5)
      XX.Fisher <- crossprod(t(XX.XW2))

      # Added by SLEE 2023/07/22
      if(!is.null(l2.var)){
        XX.Fisher[l2.var.pos, l2.var.pos] = XX.Fisher[l2.var.pos, l2.var.pos] +  1/l2.var
      }
      ##
      delta <- solve(XX.Fisher, U.star)
      #XX.covs[col.fit, col.fit] <- fast.invFisher(XX.Fisher)   ###### HERE IS THE PROBLEM!!!
    }
    if(all(is.na(XX.covs)) == T) {
      break
    }
    #delta <- as.vector(XX.covs %*% U.star)
    delta[is.na(delta)] <- 0
    mx <- max(abs(delta))/maxstep
    if (mx > 1)
      delta <- delta/mx
    evals <- evals + 1
    if (maxit > 0) {
      iter <- iter + 1
      beta <- beta + delta
      pi <- as.vector(1/(1 + exp(-x %*% beta - offset)))

    }
    if (iter == maxit | ((max(abs(delta)) <= xconv) & (all(abs(U.star[col.fit]) <
                                                           gconv))))
      break
  }
  # Error catching (if chol(x) not positive definite)
  if(all(is.na(XX.covs))==F) {

    re.out = list(beta = beta, var = var, pi = pi, hat.diag = h,
         iter = iter, evals = evals, conv = c(max(abs(U.star)),
                                              max(abs(delta))))
  }
  return(re.out)
}


# Change the algorithm, only fitting the beta[2], not intercept...
# Also no weighting (weight=1 for all)
#
# g: genotype vector
# idx: index of individuals with non-zero g
#	ex) idx = which(g>0)
fast.logistf.fit.WithL2.Sparse <- function (g, y, idx=NULL, init = NULL, out.null.firth=NULL, offset = NULL,  l2.var=NULL, control) {

  n <- length(g)
  if (is.null(offset))
    offset = rep(0, n)
  if (is.null(init))
    init = rep(0, 2)

  if (missing(control))
    control <- fast.logistf.control()

  maxit <- control$maxit
  maxstep <- control$maxstep
  maxhs <- control$maxhs
  lconv <- control$lconv
  gconv <- control$gconv
  xconv <- control$xconv

  re.out<-list(beta = NA, var = var, pi = NA, hat.diag = NA, iter = NA, evals = NA, conv = c(NA, NA, NA))

  # By SLEE, compute once...
  if(is.null(out.null.firth)){
  	out.null.firth = calculate_firth_null(y, offset=offset)
  }
  if(is.null(idx)){
  	idx=which(g>0)
  }

  pi <- out.null.firth$pi
  beta <- c(out.null.firth$beta,init[2])
  iter <- 0
  evals <- 1

  g_nozero<-g[idx]
  y_nozero<-y[idx]
  pi_nozero<-pi[idx]


  # Calculate h_nozero using the sparsity
  #
  #XW2 <- t(cbind(1,g) * ( pi * (1-pi))^0.5)
  #myQR <- qr(t(XW2))
  #Q <- qr.Q(myQR)
  #h <- (Q*Q) %*% rep(1, ncol(Q))

	XWX<-matrix(0,2,2)
 	XWX[1,1]<- out.null.firth$pi_var_sum
 	XWX[1,2]<-XWX[2,1]<-sum(g_nozero* pi_nozero * (1-pi_nozero))
 	XWX[2,2]<-sum(g_nozero^2* pi_nozero * (1-pi_nozero))

 	XW2_NoZero <- cbind(1,g_nozero) * (( pi_nozero * (1-pi_nozero))^0.5)
 	XW2_XWX_inv_NoZero =  XW2_NoZero %*% solve(XWX)
 	h_nozero <- rowSums(XW2_XWX_inv_NoZero * XW2_NoZero)




  repeat {

    beta.old <- beta
    U.star.1 <- crossprod(g_nozero,  (y_nozero - pi_nozero) + h_nozero * (0.5 - pi_nozero))

    # Added by SLEE 2023/07/22
    if(!is.null(l2.var)){
      U.star.1 = U.star.1 - 1/l2.var * beta[2]
    }

    ##
    XX.Fisher<-sum(pi_nozero * (1-pi_nozero)* g_nozero^2)

    # Added by SLEE 2023/07/22
    if(!is.null(l2.var)){
    	XX.Fisher = XX.Fisher +  1/l2.var
    }

    delta <- U.star.1 / XX.Fisher

    delta[is.na(delta)] <- 0
    mx <- max(abs(delta))/maxstep
    if (mx > 1)
      delta <- delta/mx
    evals <- evals + 1
    if (maxit > 0) {
      iter <- iter + 1
      beta[2] <- beta[2] + delta
      pi_nozero <- as.vector(1/(1 + exp(-g_nozero * beta[2] - beta[1] - offset[idx])))

    }

    if (iter == maxit | ((max(abs(delta)) <= xconv) & (all(abs(U.star.1) <
                                                           gconv))))
      break
  }


    re.out = list(beta = beta, var = var, iter = iter, evals = evals, conv = c(max(abs(U.star.1)),
                                              max(abs(delta))))

  return(re.out)
}




fast.logistf.control <- function (maxit = 200, maxhs = 15, maxstep = 15, lconv = 1e-05,
                                  gconv = 1e-05, xconv = 1e-05)
{
  list(maxit = maxit, maxhs = maxhs, maxstep = maxstep, lconv = lconv,
       gconv = gconv, xconv = xconv)
}

fast.logDet <- function (x) {
  my.chol <- tryCatch(chol(x),error=function(e) {NA})
  if (all(is.na(my.chol))==T) {
    return(NA)
  } else {
    return (2 * sum(log(diag(my.chol))))
  }
}

fast.invFisher <- function(x) {
  my.chol <- tryCatch(chol(x),error=function(e) {NA})
  if (all(is.na(my.chol))==T) {
    return(NA)
  } else {
    return (chol2inv(my.chol))
  }
  #ifelse(is.na(my.chol), NA, chol2inv(my.chol))
}


Run_Firth_With_L2<-function(G1, y, offset, l2.var, Is.Fast=TRUE, Is.Sparse=FALSE, out.null.firth=NULL, maxit=50){

  idx = NULL
  if(Is.Fast){
  	out_f1<-fast.logistf.fit.WithL2.fast(cbind(1,G1),y, offset=offset, l2.var=NULL, out.null.firth=out.null.firth)
  } else if(Is.Sparse){
  	idx = which(G1>0)
  	out_f1<-fast.logistf.fit.WithL2.Sparse(G1, y, idx=idx, offset=offset, l2.var=NULL, out.null.firth=out.null.firth)
  } else {
  	out_f1<-fast.logistf.fit.WithL2(cbind(1,G1),y, offset=offset, l2.var=NULL)
  }

  re<-out_f1$beta[2]

  # provide l2.var as regularization term
  if(!is.null(l2.var)){
  	if(Is.Fast){
    	out_f2<-fast.logistf.fit.WithL2.fast(cbind(1,G1),y, offset=offset, init=out_f1$beta, l2.var=l2.var
                                  , l2.var.pos=2, control=fast.logistf.control(maxit=maxit), out.null.firth=out.null.firth)
  	} else if(Is.Sparse){
  		out_f2<-fast.logistf.fit.WithL2.Sparse(G1, y, idx=idx, offset=offset, l2.var=l2.var, out.null.firth=out.null.firth, control=fast.logistf.control(maxit=maxit))
  	} else {
    	out_f2<-fast.logistf.fit.WithL2(cbind(1,G1),y, offset=offset, init=out_f1$beta, l2.var=l2.var
                                  , l2.var.pos=2, control=fast.logistf.control(maxit=maxit))
	}

    if(out_f2$iter <= maxit){ #Use the original version when the second one didnot converge.
      re = out_f2$beta[2]
    }
  }

  return(re)
}

Run_Firth_MultiVar_Single<-function(G1.all, obj.null, y, offset1, nMarker, l2.var=NULL, Is.Fast=TRUE, Is.Sparse=FALSE, out.null.firth=NULL){

  out_single_beta<-rep(0,nMarker)

  if(Is.Fast || Is.Sparse ){
  	if(is.null(out.null.firth)){
  		out.null.firth = calculate_firth_null(y=y, offset = offset1)
  	}
  }

  for(j in 1:nMarker){

    G1<-G1.all[,j]

    #source("Firth.R")
    re<-Run_Firth_With_L2(G1, y, offset1, l2.var, Is.Fast=Is.Fast, Is.Sparse=Is.Sparse, out.null.firth=out.null.firth)
    out_single_beta[j]<-re
  }

  # Run one more with all variants

  offset_add<-offset1 + as.vector(G1.all %*% out_single_beta)
  out_multi_beta<-rep(0,nMarker)
  for(j in 1:nMarker){

    G1<-G1.all[,j]
    offset_add1 = offset_add - G1*out_single_beta[j]

    re<-Run_Firth_With_L2(G1, y, offset_add1, l2.var,  Is.Fast=Is.Fast, Is.Sparse=Is.Sparse, out.null.firth=out.null.firth)
    out_multi_beta[j]<-re
  }

  results<-cbind(out_single_beta, out_multi_beta)
  return(results)

}


Run_Existing_Approach<-function(G1.all, obj.null, y, offset1, nMarker, l2.var=NULL){

  out_single_ext<-rep(0,nMarker)
  for(j in 1:nMarker){

    G1<-G1.all[,j]

    #source("Firth.R")
    re<-Run_Firth_With_L2(G1, y, offset1, l2.var)
    out_single_beta[j]<-re
  }

}
