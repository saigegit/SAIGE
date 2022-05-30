#Fits the null glmm
glmmkin.ai_PCG_Rcpp = function(bedFile, bimFile, famFile, Xorig, isCovariateOffset, fit0, tau=c(0,0), fixtau = c(0,0), maxiter =20, tol = 0.02, verbose = TRUE, nrun=30, tolPCG = 1e-5, maxiterPCG = 500, subPheno, indicatorGenoSamplesWithPheno, obj.noK, out.transform, tauInit, memoryChunk, LOCO, chromosomeStartIndexVec, chromosomeEndIndexVec, traceCVcutoff, isCovariateTransform, isDiagofKinSetAsOne) {
  #Fits the null generalized linear mixed model for a poisson, binomial, and gaussian
  #Args:
  #  genofile: string. Plink file for the M1 markers to be used to construct the genetic relationship matrix
  #  fit0: glm model. Logistic model output (with no sample relatedness accounted for)
  #  tau: vector for iniial values for the variance component parameter estimates
  #  fixtau: vector for fixed tau values
  #  maxiter: maximum iterations to fit the glmm model
  #  tol: tolerance for tau estimating to converge
  #  verbose: whether outputting messages in the process of model fitting
  #  nrun: integer. Number of random vectors used for trace estimation
  #  tolPCG: tolerance for PCG to converge
  #  maxiterPCG: maximum iterations for PCG to converge
  #  subPheno: data set with samples having non-missing phenotypes and non-missing genotypes (for M1 markers)
  #  obj.noK: model output from the SPAtest::ScoreTest_wSaddleApprox_NULL_Model
  #  out.transform: output from the function Covariate_Transform
  #  tauInit: vector for iniial values for the variance component parameter estimates
  #  memoryChunk: integer or float. The size (Gb) for each memory chunk
  #  LOCO:logical. Whether to apply the leave-one-chromosome-out (LOCO) option.
  #  chromosomeStartIndexVec: integer vector of length 22. Contains start indices for each chromosome, starting from 0
  #  chromosomeEndIndexVec: integer vector of length. Contains end indices for each chromosome
  #  traceCVcutoff: threshold for the coefficient of variation for trace estimation
  #Returns:
  #  model output for the null glmm

  t_begin = proc.time()
  print(t_begin)
  subSampleInGeno = subPheno$IndexGeno
  if(verbose){
    print("Start reading genotype plink file here")
  }


  re1 = system.time({setgeno(bedFile, bimFile, famFile, subSampleInGeno, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne)})
  if(verbose){
    print("Genotype reading is done")
  }

  if (LOCO){
    MsubIndVec = getQCdMarkerIndex()
    chrVec = data.table:::fread(bimFile, header = F)[,1]
    chrVec = chrVec[which(MsubIndVec == TRUE)]
    updatechrList = updateChrStartEndIndexVec(chrVec)
    LOCO = updatechrList$LOCO
    chromosomeStartIndexVec = updatechrList$chromosomeStartIndexVec
    chromosomeEndIndexVec = updatechrList$chromosomeEndIndexVec
  }

  y = fit0$y
  n = length(y)
  X = model.matrix(fit0)
  offset = fit0$offset
  if(is.null(offset)){
    offset = rep(0, n)
  }

  family = fit0$family
  eta = fit0$linear.predictors
  mu = fit0$fitted.values
  mu.eta = family$mu.eta(eta)
  Y = eta - offset + (y - mu)/mu.eta
  sqrtW = mu.eta/sqrt(fit0$family$variance(mu))
  W = sqrtW^2

  alpha0 = fit0$coef
  eta0 = eta

  if(family$family %in% c("poisson", "binomial")) {
    tau[1] = 1
    fixtau[1] = 1


    if(tauInit[fixtau == 0] == 0){
      tau[fixtau == 0] = 0.1
    }else{
      tau[fixtau == 0] = tauInit[fixtau == 0]
    }
    q = 1
    cat("inital tau is ", tau,"\n")
    tau0=tau
    re.coef = Get_Coef(y, X, tau, family, alpha0, eta0,  offset,verbose=verbose, maxiterPCG=maxiterPCG, tolPCG = tolPCG, maxiter=maxiter)
    re = getAIScore(re.coef$Y, X, re.coef$W, tau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov, nrun, maxiterPCG,tolPCG = tolPCG, traceCVcutoff = traceCVcutoff)
    tau[2] = max(0, tau0[2] + tau0[2]^2 * (re$YPAPY - re$Trace)/n)

  }else if(family$family == "gaussian"){
    if(sum(tauInit[fixtau == 0]) == 0){
      #tau[fixtau == 0] = var(Y)/(q+1)
      tau[1] = 1
      tau[2] = 0
      if (abs(var(Y)) < 0.1){
        stop("WARNING: variance of the phenotype is much smaller than 1. Please consider invNormalize=T\n")
      }
    }else{
      tau[fixtau == 0] = tauInit[fixtau == 0]
    }
    q = 1
    cat("inital tau is ", tau,"\n")
    tau0=tau
    re.coef = Get_Coef(y, X, tau, family, alpha0, eta0,  offset,verbose=verbose, maxiterPCG=maxiterPCG, tolPCG = tolPCG, maxiter=maxiter) 
    re = getAIScore_q(re.coef$Y, X, re.coef$W, tau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov, nrun, maxiterPCG,tolPCG = tolPCG, traceCVcutoff = traceCVcutoff)
    tau[2] = max(0, tau0[2] + tau0[2]^2 * (re$YPAPY - re$Trace[2])/n)
    tau[1] = max(0, tau0[1] + tau0[1]^2 * (re$YPA0PY - re$Trace[1])/n)
  }


  if(verbose) {
    cat("Variance component estimates:\n")
    print(tau)
  }

  for (i in seq_len(maxiter)) {
    W = sqrtW^2

    if(verbose) cat("\nIteration ", i, tau, ":\n")
      alpha0 = re.coef$alpha
      tau0 = tau
      cat("tau0_v1: ", tau0, "\n")
      eta0 = eta
      # use Get_Coef before getAIScore
      t_begin_Get_Coef = proc.time()
      re.coef = Get_Coef(y, X, tau, family, alpha0, eta0,  offset, verbose=verbose, maxiterPCG=maxiterPCG, tolPCG = tolPCG, maxiter=maxiter)
      t_end_Get_Coef =  proc.time()
      cat("t_end_Get_Coef - t_begin_Get_Coef\n")
      print(t_end_Get_Coef - t_begin_Get_Coef)
      fit = fitglmmaiRPCG(re.coef$Y, X, re.coef$W, tau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov, nrun, maxiterPCG, tolPCG, tol = tol, traceCVcutoff = traceCVcutoff)
      t_end_fitglmmaiRPCG= proc.time()
      cat("t_end_fitglmmaiRPCG - t_begin_fitglmmaiRPCG\n")
      print(t_end_fitglmmaiRPCG - t_end_Get_Coef)

      tau = as.numeric(fit$tau)
      cov = re.coef$cov
      alpha = re.coef$alpha
      eta = re.coef$eta
      Y = re.coef$Y
      mu = re.coef$mu

     
      print(abs(tau - tau0)/(abs(tau) + abs(tau0) + tol))
      cat("tau: ", tau, "\n")
      cat("tau0: ", tau0, "\n")

      if(family$family == "gaussian"){
        if(tau[1]<=0){
          stop("ERROR! The first variance component parameter estimate is 0\n")
        }
      }

      if(tau[2] == 0) break
      
      # Use only tau for convergence evaluation, because alpha was evaluated already in Get_Coef
      if(max(abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol) break

      if(max(tau) > tol^(-2)) {
        warning("Large variance estimate observed in the iterations, model not converged...", call. = FALSE)
        i = maxiter
        break
      }
  }

  if(verbose) cat("\nFinal " ,tau, ":\n")

    #added these steps after tau is estimated 04-14-2018
  re.coef = Get_Coef(y, X, tau, family, alpha, eta,  offset,verbose=verbose, maxiterPCG=maxiterPCG, tolPCG = tolPCG, maxiter=maxiter)
  cov = re.coef$cov
  alpha = re.coef$alpha
  eta = re.coef$eta
  Y = re.coef$Y
  mu = re.coef$mu

  converged = ifelse(i < maxiter, TRUE, FALSE)
  res = y - mu

  if(family$family == "binomial"){
    mu2 = mu * (1-mu)
    traitType = "binary"
  }else if(family$family == "poisson"){
    mu2 = mu
    traitType = "count"
  }else if(family$family == "gaussian"){
    mu2 = rep((1/(tau[1])),length(res))
    traitType = "quantitative"
  }

  #if(isCovariateTransform & hasCovariate){
  if(!is.null(out.transform) & is.null(fit0$offset)){
    coef.alpha<-Covariate_Transform_Back(alpha, out.transform$Param.transform)
  }else{
    coef.alpha = alpha
  }

  #mu2 = mu * (1-mu)

  if(!isCovariateOffset){
    obj.noK = ScoreTest_NULL_Model(mu, mu2, y, X)
    glmmResult = list(theta=tau, coefficients=coef.alpha, linear.predictors=eta, fitted.values=mu, Y=Y, residuals=res, cov=cov, converged=converged,sampleID = subPheno$IID, obj.noK=obj.noK, y = y, X = X, traitType=traitType, isCovariateOffset = isCovariateOffset)
  }else{
    obj.noK = ScoreTest_NULL_Model(mu, mu2, y, Xorig)
    glmmResult = list(theta=tau, coefficients=coef.alpha, linear.predictors=eta, fitted.values=mu, Y=Y, residuals=res, cov=cov, converged=converged,sampleID = subPheno$IID, obj.noK=obj.noK, y = y, X = Xorig, traitType=traitType, isCovariateOffset = isCovariateOffset)
  }


  #LOCO: estimate fixed effect coefficients, random effects, and residuals for each chromoosme

  glmmResult$LOCO = LOCO
  t_end_null = proc.time()
  cat("t_end_null - t_begin, fitting the NULL model without LOCO took\n")
  print(t_end_null - t_begin)
  if(LOCO){
    set_Diagof_StdGeno_LOCO()
    glmmResult$LOCOResult = list()
    for (j in 1:22){
      startIndex = chromosomeStartIndexVec[j]
      endIndex = chromosomeEndIndexVec[j]
      if(!is.na(startIndex) && !is.na(endIndex)){
        cat("leave chromosome ", j, " out\n")
        setStartEndIndex(startIndex, endIndex, j-1)
        t_begin_Get_Coef_LOCO = proc.time()
        re.coef_LOCO = Get_Coef_LOCO(y, X, tau, family, alpha, eta,  offset,verbose=verbose, maxiterPCG=maxiterPCG, tolPCG = tolPCG, maxiter=maxiter)
        t_end_Get_Coef_LOCO = proc.time()
        cat("t_end_Get_Coef_LOCO - t_begin_Get_Coef_LOCO\n")
        print(t_end_Get_Coef_LOCO - t_begin_Get_Coef_LOCO)
        cov = re.coef_LOCO$cov
        alpha = re.coef_LOCO$alpha
        eta = re.coef_LOCO$eta
        Y = re.coef_LOCO$Y
        mu = re.coef_LOCO$mu
        #mu2 = mu * (1-mu)
        #mu2 = mu
        res = y - mu

        if(family$family == "binomial"){
          mu2 = mu * (1-mu)
        }else if(family$family == "poisson"){
          mu2 = mu
        }else if(family$family == "gaussian"){
          mu2 = rep((1/(tau[1])),length(res))
        }

        if(!is.null(out.transform) & is.null(fit0$offset)){
          coef.alpha<-Covariate_Transform_Back(alpha, out.transform$Param.transform)
        }else{
          coef.alpha = alpha
        }
        if(!isCovariateOffset){
          obj.noK = ScoreTest_NULL_Model(mu, mu2, y, X)
        }else{
          obj.noK = ScoreTest_NULL_Model(mu, mu2, y, Xorig)
        }
        glmmResult$LOCOResult[[j]] = list(isLOCO = TRUE, coefficients=coef.alpha, linear.predictors=eta, fitted.values=mu, Y=Y, residuals=res, cov=cov, obj.noK = obj.noK)
      }else{
        glmmResult$LOCOResult[[j]] = list(isLOCO = FALSE)
      }
    }
  }

  return(glmmResult)
}
