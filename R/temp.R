  family = obj.glm.null$family
    print(family)
      eta = obj.glmm.null$linear.predictors
        mu = obj.glmm.null$fitted.values
	  mu.eta = family$mu.eta(eta)
	    sqrtW = mu.eta/sqrt(obj.glm.null$family$variance(mu))
	      W = sqrtW^2   ##(mu*(1-mu) for binary)
	        tauVecNew = obj.glmm.null$theta
		  X = obj.glmm.null$X


		    Sigma_iX_noLOCO = getSigma_X(W, tauVecNew, X, maxiterPCG, tolPCG)



		    1124



		     glmmResult = list(theta=tau, coefficients=coef.alpha, linear.predictors=eta, fitted.values=mu, Y=Y, residuals=res, cov=cov, converged=converged,sampleID = subPheno$IID, obj.noK=obj.noK, y = y, X = Xorig, traitType="binary", isCovariateOffset = isCovariateOffset)


glmmkin.ai_PCG_Rcpp_Quantitative



useSparseGRMtoFitNULL


m_Sigma_iXXSigma_iX
