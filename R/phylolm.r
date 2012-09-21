phylolm <- function(formula, data=list(), phy, 
	model=c("BM","OUrandomRoot","OUfixedRoot","lambda","kappa","delta","EB","trend"),
	lower.bound=NULL, upper.bound=NULL, starting.value=NULL, ...) 
{
	#----------
	### initialize	
	if (!inherits(phy, "phylo")) stop("object \"phy\" is not of class \"phylo\".")
	model = match.arg(model)	
	if ((model=="trend")&(is.ultrametric(phy)))
		stop("can not estimate trend for ultrametric tree.")	
      if (is.null(phy$edge.length)) stop("the tree has no branch lengths.")
	if (is.null(phy$tip.label)) stop("the tree has no tip labels.")	
	tol = 1e-10	
	phy = reorder(phy,"pruningwise")
	n <- length(phy$tip.label)
    	N <- dim(phy$edge)[1]
    	ROOT <- n + 1L
	anc <- phy$edge[, 1]
    	des <- phy$edge[, 2]
	mf = model.frame(formula=formula,data=data)
	if (nrow(mf)!=length(phy$tip.label)) stop("number of rows of the design matrix does not match with length of the tree.")
	if (is.null(rownames(mf))) {
		cat("Warning: no tip labels, order assumed to be the same as in the tree.\n")
		data.names = phy$tip.label 
		}
	else data.names = rownames(mf)
	order = match(data.names, phy$tip.label)
	if (sum(is.na(order))>0) {
		cat("Warning: data names do not match with the tip labels.\n")
		rownames(mf) = data.names
		} else {
			tmp = mf
			rownames(mf) = phy$tip.label
			mf[order,] = tmp[1:nrow(tmp),]
		}
	X = model.matrix(attr(mf, "terms"), data=mf)
 	y = model.response(mf)
	d = ncol(X)
	OU = c("OUrandomRoot","OUfixedRoot")
	flag = 0

	### Default bounds
	bounds.default = matrix(c(1e-7,50,1e-7,1,1e-6,1,1e-5,3,-3,0), ncol=2, byrow=TRUE)
      rownames(bounds.default) = c("alpha","lambda","kappa","delta","rate")
      colnames(bounds.default) = c("min","max")

	### Default starting values
	starting.values.default = c(0.5,0.5,0.5,0.5,-1)
      names(starting.values.default) = c("alpha","lambda","kappa","delta","rate")

	### User defined bounds and starting values
	if (is.null(lower.bound)) {
		if (model %in% OU) 
			lower.bound = bounds.default[1,1]
		if (model=="lambda") lower.bound = bounds.default[2,1]
		if (model=="kappa") lower.bound = bounds.default[3,1]
		if (model=="delta") lower.bound = bounds.default[4,1]
		if (model=="EB") lower.bound = bounds.default[5,1]
	}
	if (is.null(upper.bound)) {
		if (model %in% OU) 
			upper.bound = bounds.default[1,2]
		if (model=="lambda") upper.bound = bounds.default[2,2]
		if (model=="kappa") upper.bound = bounds.default[3,2]
		if (model=="delta") upper.bound = bounds.default[4,2]
		if (model=="EB") upper.bound = bounds.default[5,2]
	}
	if (is.null(starting.value)) {
		if (model %in% OU) 
			starting.value = starting.values.default[1]
		if (model=="lambda") starting.value = starting.values.default[2]
		if (model=="kappa") starting.value = starting.values.default[3]
		if (model=="delta") starting.value = starting.values.default[4]
		if (model=="EB") starting.value = starting.values.default[5]
	}
	#----------

	#----------
	### branching.time function for ultrametric tree in pruningwise order
	pruningwise.branching.times <- function(phy) {   	
   	xx <- numeric(phy$Nnode)
      interns <- which(phy$edge[, 2] > n)
    	for (j in length(interns):1) {
		i = interns[j]
		xx[des[i] - n] <- xx[anc[i] - n] + phy$edge.length[i]
    		}
    	depth <- xx[phy$edge[1, 1] - n] + phy$edge.length[1]
    	xx <- depth - xx
    	names(xx) <- if (is.null(phy$node.label)) (n + 1):(n + phy$Nnode) else phy$node.label
      return(xx)
      }
	#----------

	#----------
	### distance from the root of tree's nodes in pruningwise order
	pruningwise.distFromRoot <- function(phy) {   
	xx <- numeric(phy$Nnode+n)
      for (i in N:1) xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]
      names(xx) <- if (is.null(phy$node.label)) 1:(n + phy$Nnode) else phy$node.label
      return(xx)
 	}
	#----------

	#----------
	### preparing for OU model
	if (model %in% OU) {
		if (!is.ultrametric(phy)) flag = 1
		dist = pruningwise.distFromRoot(phy)
		D = max(dist[1:n]) - dist[1:n]
		D = D - mean(D)
		for (i in 1:n) {
			tmp = phy$edge.length[which(des==i)]
			phy$edge.length[which(des==i)] = tmp + D[i]			
			}
		times <- pruningwise.branching.times(phy)
		names(times) <- (n+1):(n+phy$Nnode) 
  		Tmax <- max(times)
		}
	#----------
	
	#----------
	### preparing for trend model
	if (model == "trend") {
		trend = pruningwise.distFromRoot(phy)[1:n]
		X = cbind(X,trend)
		d = d+1
		}	
	#----------

	#----------
	### function to calculate transformed branch lengths
	transformed.branch.lengths <- function(model, value = NULL)
	{
		### BM model
		if (model %in% c("BM","trend")) {
			edge.length = phy$edge.length
			root.edge = 0
			}
		
		### OUrandomRoot model	
		if (model=="OUrandomRoot") {
			alpha = value
			edge.length = numeric(N)
			distFromRoot <-  exp(-2*alpha*times)
			externalEdge <- des<=n
			for (i in 1:N) {
				d1 = distFromRoot[which(names(times)==anc[i])]
				if (externalEdge[i]) {d2 = exp(-2*alpha*D[des[i]])} else d2 = distFromRoot[which(names(times)==des[i])]
				edge.length[i] = d2 - d1
				}
			root.edge = min(distFromRoot)
			}

		### OUfixedRoot model	
		if (model=="OUfixedRoot") {
			alpha = value
			edge.length = numeric(N)
			distFromRoot <-  exp(-2*alpha*times)*(1 - exp(-2*alpha*(Tmax-times)))
			externalEdge <- des<=n
			for (i in 1:N) {
				d1 = distFromRoot[which(names(times)  == anc[i])]
				if (externalEdge[i]) {d2 = exp(-2*alpha*D[des[i]])*(1-exp(-2*alpha*(Tmax-D[des[i]])))}
					else d2 = distFromRoot[which(names(times)==des[i])]
				edge.length[i] = d2 - d1
				}
			root.edge = min(distFromRoot)
			}
	
		### lambda model
		if (model=="lambda") {
			lambda = value
			edge.length = numeric(N)
			distFromRoot <- pruningwise.distFromRoot(phy)
			externalEdge <- des<=n
			for (i in 1:N) {				
				if (externalEdge[i]) {edge.length[i]  = phy$edge.length[i]*lambda+(1-lambda)*distFromRoot[des[i]]}
					else edge.length[i] = phy$edge.length[i]*lambda
				}
			root.edge = 0
			}
		### kappa model
		if (model=="kappa") {
			kappa = value
			edge.length = phy$edge.length^kappa
			root.edge = 0
			}

		### delta model
		if (model=="delta") {
			delta = value
			edge.length = numeric(N)
			distFromRoot <- pruningwise.distFromRoot(phy)
			depth = max(distFromRoot)
			for (i in 1:N) edge.length[i] = (distFromRoot[des[i]]^delta - distFromRoot[anc[i]]^delta)*depth^(1-delta)
			root.edge = 0
			}

		### EB model
		if (model=="EB") {
			rate = value
			edge.length = numeric(N)
			if (rate==0) edge.length = phy$edge.length
				else {
					distFromRoot <- pruningwise.distFromRoot(phy)
					for (i in 1:N) edge.length[i] = (exp(rate*distFromRoot[des[i]])-exp(rate*distFromRoot[anc[i]]))/rate
					}			
			root.edge = 0
			}
	return(list(edge.length,root.edge))
	}
	#----------

	#----------	
	### a function using BM structure to do computation
	BMstructure.compute <- function(trans,y,X) {		
	zero <- rep(0,n + phy$Nnode)
	vec11 <- rep(0,n + phy$Nnode) # 1'V^{-1}1
	y1 <- rep(0,n + phy$Nnode) # y'V^{-1}1
	yy <- rep(0,n + phy$Nnode) # y'V^{-1}y
	logd <- rep(0,n + phy$Nnode) # logdetV
	X1 <- matrix(0,d,n + phy$Nnode) # X'V^{-1}1
	XX <- array(0,dim=c(d,d,n + phy$Nnode)) # X'V^{-1}X
	Xy <- matrix(0,d,n + phy$Nnode) # X'V^{-1}y
	el = trans[[1]]
	root.edge = trans[[2]]
	for (i in 1:N) {
		if (des[i] <= n) {		
			if (el[i] > 0) {
		   		logd[des[i]] = log(el[i])
		   		Xy[,des[i]] = X[des[i],]*y[des[i]]/el[i]
		   		yy[des[i]] = y[des[i]]^2/el[i]
		   		y1[des[i]] = y[des[i]]/el[i]
		   		XX[,,des[i]] = X[des[i],]%*%t(X[des[i],])/el[i]
		   		X1[,des[i]] = X[des[i],]/el[i]
		   		vec11[des[i]] = 1/el[i] 		  
				} else zero[anc[i]] = zero[anc[i]]+1				
			} else {
			if (((el[i]<=0) & (zero[des[i]]>0)) | (zero[des[i]]>1))
				stop("logLikelihood = -Inf.")	
			logd[des[i]] = logd[des[i]] + log(1+el[i]*vec11[des[i]])
			Xy[,des[i]] = Xy[,des[i]] - el[i]/(1+el[i]*vec11[des[i]])*X1[,des[i]]*y1[des[i]]
		   	yy[des[i]] = yy[des[i]] - el[i]/(1+el[i]*vec11[des[i]])*y1[des[i]]^2
		   	y1[des[i]] = y1[des[i]]*(1 - 1/(1+1/el[i]/vec11[des[i]]))
			XX[,,des[i]] = XX[,,des[i]] - el[i]/(1+el[i]*vec11[des[i]])*X1[,des[i]]%*%t(X1[,des[i]])
		   	X1[,des[i]] = X1[,des[i]]*(1 - 1/(1+1/el[i]/vec11[des[i]]))		   	
			vec11[des[i]] = vec11[des[i]]*(1 - 1/(1+1/el[i]/vec11[des[i]]))
		   	}
   		logd[anc[i]] = logd[anc[i]] + logd[des[i]]
		Xy[,anc[i]] = Xy[,anc[i]] + Xy[,des[i]]
		yy[anc[i]] = yy[anc[i]] + yy[des[i]]
		y1[anc[i]] = y1[anc[i]] + y1[des[i]]
		XX[,,anc[i]] = XX[,,anc[i]] + XX[,,des[i]]
		X1[,anc[i]] = X1[,anc[i]] + X1[,des[i]]
		vec11[anc[i]] = vec11[anc[i]] + vec11[des[i]]				
	}
	logd[ROOT] = logd[ROOT] + log(1+root.edge*vec11[ROOT])
	Xy[,ROOT] = Xy[,ROOT] - root.edge/(1+root.edge*vec11[ROOT])*X1[,ROOT]*y1[ROOT]
	yy[ROOT] = yy[ROOT] - root.edge/(1+root.edge*vec11[ROOT])*y1[ROOT]^2
	y1[ROOT] = y1[ROOT]*(1 - 1/(1+1/root.edge/vec11[ROOT]))
	XX[,,ROOT] = XX[,,ROOT] - root.edge/(1+root.edge*vec11[ROOT])*X1[,ROOT]%*%t(X1[,ROOT])
	X1[,ROOT] = X1[,ROOT]*(1 - 1/(1+1/root.edge/vec11[ROOT]))
	vec11[ROOT] = vec11[ROOT]*(1 - 1/(1+1/root.edge/vec11[ROOT]))
	comp = list(vec11=vec11[ROOT],y1=y1[ROOT],yy=yy[ROOT],X1=X1[,ROOT],XX=XX[,,ROOT],Xy=Xy[,ROOT],logd=logd[ROOT]) 
	return(comp)
	}
	#----------
	
	#----------
	### Estimate with BM stucture
	estimate.BMstructure <- function(trans,y,X) {
	comp = BMstructure.compute(trans,y,X)
	invXX = solve(comp$XX)	
	betahat = invXX%*%comp$Xy
	sigma2hat = as.numeric((comp$yy - 2*t(betahat)%*%comp$Xy + t(betahat)%*%comp$XX%*%betahat)/n)
	if (sigma2hat<0) stop("The estimate of sigma square is negative.")
	llh = as.numeric(- n/2*log(2*pi) - n/2 - n/2*log(sigma2hat) - comp$logd/2)
	vcov = sigma2hat*invXX*n/(n-d)
	return(list(llh=llh, betahat = as.vector(betahat), sigma2hat=sigma2hat,vcov=vcov))
	}
	#----------

	#----------
	### Optimization routine
	opt.routine <- function(model,lower,upper,start) {
	if ((lower>start)||(upper<start))
		stop("The starting value is not within the bounds of the parameter.")
	minus.pfllh <- function(model, logvalue) {
		if (model == "EB") value=logvalue else value=exp(logvalue)
		trans = transformed.branch.lengths(model, value = value)
		if (flag == 1) {
			BMest = estimate.BMstructure(trans, y=exp(-value*D)*y, X=exp(-value*D)*X)
			return(-BMest$llh + value*sum(D)) 
			} else {
				BMest = estimate.BMstructure(trans, y=y, X=X)
 				return(-BMest$llh) 
				}
		}
	if (lower==upper) {
		MLEvalue = lower
		} else {
			if (model == "EB") {
				logstart = start
				loglower = lower
				logupper = upper 
				} else {
					logstart = log(start)
					loglower = log(lower)
					logupper = log(upper)
					}
			opt <- optim(par = logstart, 
              		fn = function(par){minus.pfllh(model, logvalue = par)},
              		method = "L-BFGS-B",lower=loglower, upper = logupper, ...)
			if (model == "EB") MLEvalue = as.numeric(opt$par) else MLEvalue = as.numeric(exp(opt$par))
			}
	if ((isTRUE(all.equal(MLEvalue,lower, tol=tol)))||(isTRUE(all.equal(MLEvalue,upper,tol=tol)))) {
		if (model %in% OU) parameter = "alpha"
		if (model %in% c("lambda","kappa","delta")) parameter = model
		if (model == "EB") parameter = "rate"
		matchbound = 1
		if ((model %in% c("lambda","kappa"))&&(MLEvalue == 1)) matchbound=0
		if ((model == "EB")&&(MLEvalue == 0)) matchbound=0
		if (matchbound == 1) cat("Warning: the estimation of", parameter, "matches the upper/lower bound for this parameter.\n You may change the bounds using options \"upper.bound\" and \"lower.bound\".\n")
		}
	trans = transformed.branch.lengths(model, value = MLEvalue)
	if (flag == 1) {
		BMest = estimate.BMstructure(trans, y=exp(-MLEvalue*D)*y, X=exp(-MLEvalue*D)*X)
		return(list(llh=BMest$llh-MLEvalue*sum(D), betahat=BMest$betahat, sigma2hat=BMest$sigma2hat, MLEvalue=MLEvalue, vcov = BMest$vcov))
		} else {
			BMest = estimate.BMstructure(trans, y=y, X=X)
			return(list(llh=BMest$llh, betahat=BMest$betahat, sigma2hat=BMest$sigma2hat, MLEvalue=MLEvalue, vcov = BMest$vcov))
			}	
	}
	#----------

	#------------------------------
	### Fitting
	lower = lower.bound
	upper = upper.bound
	start = starting.value

	if (model %in% c("BM","trend")) {
		trans = transformed.branch.lengths(model)	
		BMest = estimate.BMstructure(trans, y=y, X=X)
		results = list(coefficients=BMest$betahat, sigma2=BMest$sigma2hat, 
			optpar=NULL, logLik=BMest$llh,
			p=1+d, aic=2*(1+d)-2*BMest$llh, vcov = BMest$vcov)
		}

	if (model %in% OU) {	
		opt = opt.routine(model,lower,upper,start)
		results = list(coefficients=opt$betahat, sigma2=2*opt$MLEvalue*opt$sigma2hat, 
			optpar=opt$MLEvalue, logLik=opt$llh,
			p=2+d, aic=2*(2+d)-2*opt$llh, vcov = opt$vcov)
		}

	if (model=="lambda") {
		opt = opt.routine(model,lower,upper,start)
		results = list(coefficients=opt$betahat, sigma2=opt$sigma2hat, 
			optpar=opt$MLEvalue, logLik=opt$llh, 
			p=2+d, aic=2*(2+d)-2*opt$llh, vcov = opt$vcov)		
		}

	if (model=="kappa") {
		opt = opt.routine(model,lower,upper,start)
		results = list(coefficients=opt$betahat, sigma2=opt$sigma2hat, 
			optpar=opt$MLEvalue, logLik=opt$llh, 
			p=2+d, aic=2*(2+d)-2*opt$llh, vcov = opt$vcov)
		}

	if (model=="delta") {
		opt = opt.routine(model,lower,upper,start)
		results = list(coefficients=opt$betahat, sigma2=opt$sigma2hat, 
			optpar=opt$MLEvalue, logLik=opt$llh, 
			p=2+d, aic=2*(2+d)-2*opt$llh, vcov = opt$vcov)
		}

	if (model=="EB") {	
		opt = opt.routine(model,lower,upper,start)
		results = list(coefficients=opt$betahat, sigma2=opt$sigma2hat, 
			optpar=opt$MLEvalue, logLik=opt$llh, 
			p=2+d, aic=2*(2+d)-2*opt$llh, vcov = opt$vcov)
		}

	#------------------------------

	names(results$coefficients) = colnames(X)
	colnames(results$vcov) = colnames(X)
	rownames(results$vcov) = colnames(X)
	results$fitted.values = X %*% results$coefficients
	results$residuals = y - results$fitted.values
	results$y = y
	results$n = n
	results$d = d
	results$formula = formula
	results$call = match.call()
	results$model = model
	class(results) = "phylolm"
	results
}

print.phylolm <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
	cat("Call:\n")
	print(x$call)
	cat("\n")
	aiclogLik = c(x$aic,x$logLik)
	names(aiclogLik) = c("AIC","logLik")
	print(aiclogLik, digits = digits)
	cat("\nParameter estimate(s) using ML method:\n")
	if (!is.null(x$optpar)) {
		if (x$model %in% c("OUrandomRoot","OUfixedRoot")) cat("alpha:",x$optpar)
		if (x$model %in% c("lambda","kappa","delta")) cat(x$model,":",x$optpar)
		if (x$model=="EB") cat("rate:",x$optpar)
		cat("\n")
		}
	cat("sigma2:",x$sigma2)
	cat("\n")

	cat("\nCoefficients:\n")
	print(x$coefficients)
}

summary.phylolm <- function(object, ...) 
{
	se <- sqrt(diag(object$vcov))
	tval <- coef(object) / se
	TAB <- cbind(Estimate = coef(object),
		StdErr = se,
		t.value = tval,
		p.value = 2*pt(-abs(tval), df=object$n - object$d))
	res <- list(call=object$call,
		coefficients=TAB,
		residuals = object$residuals,
		sigma2 = object$sigma2,
		optpar=object$optpar,
		logLik=object$logLik,
		df=object$p,
		aic=object$aic,
		model=object$model)
	class(res) = "summary.phylolm"
	res
}

print.summary.phylolm <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
	cat("\nCall:\n")
	print(x$call)
	cat("\n")
	aiclogLik = c(x$aic,x$logLik)
	names(aiclogLik) = c("AIC","logLik")
	print(aiclogLik, digits = digits)
	r <- zapsmall(quantile(x$residuals), digits + 1)
	names(r) <- c("Min", "1Q", "Median", "3Q", "Max")
	cat("\nRaw residuals:\n")
	print(r, digits = digits)

	cat("\nParameter estimate(s) using ML method:\n")
	if (!is.null(x$optpar)) {
		if (x$model %in% c("OUrandomRoot","OUfixedRoot")) cat("alpha:",x$optpar)
		if (x$model %in% c("lambda","kappa","delta")) cat(x$model,":",x$optpar)
		if (x$model=="EB") cat("rate:",x$optpar)
		cat("\n")
		}
	cat("sigma2:",x$sigma2)
	cat("\n")
	cat("\nCoefficients:\n")
	printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
	if (!is.null(x$optpar)) {
		cat("\nNote: p-values are conditional on ")
		if (x$model %in% c("OUrandomRoot","OUfixedRoot")) cat("alpha=",x$optpar,".",sep="")
		if (x$model %in% c("lambda","kappa","delta")) cat(x$model,"=",x$optpar,".",sep="")
		if (x$model=="EB") cat("rate=",x$optpar,".",sep="")
		}
	cat("\n")
}

residuals.phylolm <-function(object,type=c("response"), ...)
{
	type <- match.arg(type)
    	r <- object$residuals
	r	 
}

vcov.phylolm <- function(object, ...)
{
	vcov = object$vcov
	vcov
}

logLik.phylolm <- function(object, ...)
{
	res = list(logLik = object$logLik, df = object$p)
	class(res) = "logLik.phylolm"
	res
}

print.logLik.phylolm <- function (x, ...) 
{
	cat("'log Lik.' ",x$logLik," (df=",x$df,")\n", sep = "")
}

predict.phylolm <- function(object, newdata=NULL, ...)
{
	if (object$model=="trend")
		stop("Predicting for trend model has not been implemented.")
	if(is.null(newdata)) y <- fitted(object)
		else{			
			X = model.matrix(delete.response(terms(formula(object))),data = newdata)
			y <- X %*% coef(object)
		}
	y
}

plot.phylolm <-function(x, ...)
{
	plot(x$y, fitted(x), xlab = "Observed value", ylab = "Fitted value")
}

rTrait <- function(phy, model=c("BM","OU","lambda","kappa","delta","EB","trend"),
	parameters=NULL)
{
	#----------
	### initialize
	if (!inherits(phy, "phylo")) stop("object \"phy\" is not of class \"phylo\".")
	model = match.arg(model)
      if (is.null(phy$edge.length)) stop("the tree has no branch lengths.")
	if (is.null(phy$tip.label)) stop("the tree has no tip labels.")	
	if ((model=="trend")&(is.ultrametric(phy)))
		stop("should not use trend model for ultrametric tree.")
	phy = reorder(phy,"pruningwise")
	n <- length(phy$tip.label)
    	N <- dim(phy$edge)[1]
    	ROOT <- n + 1L
	anc <- phy$edge[, 1]
    	des <- phy$edge[, 2]

	### Default parameters
	parameters.default = c(0,1,0,0,1,1,1,0,0)
	names(parameters.default) = c("ancestral.state", "sigma2", 
		"optimal.value", "alpha", "lambda", "kappa", "delta", "rate", "trend")

	### User defined parameters
	if (is.null(parameters)) {
		parameters = parameters.default
		} else { 
			if (class(parameters)!= "list") {
				stop("please specify parameters as a list().")
				} else {
					specified <- !c(is.null(parameters$ancestral.state), 
						is.null(parameters$sigma2), 
						is.null(parameters$optimal.value),
						is.null(parameters$alpha),
						is.null(parameters$lambda),
						is.null(parameters$kappa),
						is.null(parameters$delta),
						is.null(parameters$rate),
						is.null(parameters$trend))
					parameters.user <- c(parameters$ancestral.state,
						parameters$sigma2,
						parameters$optimal.value,
						parameters$alpha,
						parameters$lambda,
						parameters$kappa,
						parameters$delta,
						parameters$rate,
						parameters$trend)
					names(parameters.default) = c("ancestral.state","sigma2","optimal.value",
						"alpha", "lambda", "kappa", "delta", "rate", "trend")
					parameters <- parameters.default
		 			parameters[specified] <- parameters.user 
					}				
			}
	p = list(ancestral.state = parameters[1],
		sigma2 = parameters[2],
		optimal.value = parameters[3],
		alpha = parameters[4],
		lambda = parameters[5],
		kappa = parameters[6],
		delta = parameters[7],
		rate = parameters[8],
		trend = parameters[9])
	#----------
	
	#----------
	### branching.time function for ultrametric tree in pruningwise order
	pruningwise.branching.times <- function(phy) {   	
   	xx <- numeric(phy$Nnode)
      interns <- which(phy$edge[, 2] > n)
    	for (j in length(interns):1) {
		i = interns[j]
		xx[des[i] - n] <- xx[anc[i] - n] + phy$edge.length[i]
    		}
    	depth <- xx[phy$edge[1, 1] - n] + phy$edge.length[1]
    	xx <- depth - xx
    	names(xx) <- if (is.null(phy$node.label)) (n + 1):(n + phy$Nnode) else phy$node.label
      return(xx)
      }
	#----------

	#----------
	### distance from the root of tree's nodes in pruningwise order
	pruningwise.distFromRoot <- function(phy) {   
	xx <- numeric(phy$Nnode+n)
      for (i in N:1) xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]
      names(xx) <- if (is.null(phy$node.label)) 1:(n + phy$Nnode) else phy$node.label
      return(xx)
 	}
	#----------

	#----------
	### function to calculate transformed branch lengths
	transformed.branch.lengths <- function(model, value=NULL)
	{
		### BM model
		if (model %in% c("BM","trend")) {
			edge.length = phy$edge.length
			root.edge = 0
			}
		
		### lambda model
		if (model=="lambda") {
			lambda = value
			edge.length = numeric(N)
			distFromRoot <- pruningwise.distFromRoot(phy)
			externalEdge <- des<=n
			for (i in 1:N) {				
				if (externalEdge[i]) {edge.length[i]  = phy$edge.length[i]*lambda+(1-lambda)*distFromRoot[des[i]]}
					else edge.length[i] = phy$edge.length[i]*lambda
				}
			root.edge = 0
			}
		### kappa model
		if (model=="kappa") {
			kappa = value
			edge.length = phy$edge.length^kappa
			root.edge = 0
			}

		### delta model
		if (model=="delta") {
			delta = value
			edge.length = numeric(N)
			distFromRoot <- pruningwise.distFromRoot(phy)
			depth = max(distFromRoot)
			for (i in 1:N) edge.length[i] = (distFromRoot[des[i]]^delta - distFromRoot[anc[i]]^delta)*depth^(1-delta)
			root.edge = 0
			}

		### EB model
		if (model=="EB") {
			rate = value
			edge.length = numeric(N)
			if (rate==0) edge.length = phy$edge.length
				else {
					distFromRoot <- pruningwise.distFromRoot(phy)
					for (i in 1:N) edge.length[i] = (exp(rate*distFromRoot[des[i]])-exp(rate*distFromRoot[anc[i]]))/rate
					}			
			root.edge = 0
			}

	return(list(edge.length,root.edge))
	}
	#----------
	
	#----------
	### Simulation from BM model
	BMsim <- function(trans, ancestral.state, sigma2) {
		x <- numeric(n + phy$Nnode)	
		el = trans[[1]]
		root.edge = trans[[2]]
      	x[ROOT] <- rnorm(1,ancestral.state,sqrt(root.edge*sigma2))       
      	for (i in N:1) x[des[i]] = rnorm(1,x[anc[i]],sqrt(el[i]*sigma2))
      	x <- x[1:n]
      	names(x) <- phy$tip.label
      	return(x)
	}
	### Simulation from OU model
	OUsim <- function(y0, mu, alpha, gam) {
		el = phy$edge.length
		x <- numeric(n + phy$Nnode)	
	  	x[ROOT] = y0
	      for (i in N:1) x[des[i]] = rnorm(1,x[anc[i]]*exp(-alpha*el[i])+mu*(1-exp(-alpha*el[i])),sqrt(gam*(1-exp(-2*alpha*el[i]))))
	      x <- x[1:n]
		names(x) <- phy$tip.label
		return(x)
	}
	### Simulation from BM model with a trend
	BMtrend.sim <- function (trans, ancestral.state, sigma2, trend) {
		x <- numeric(n + phy$Nnode)	
		el = trans[[1]]
		root.edge = trans[[2]]
      	x[ROOT] <- rnorm(1,ancestral.state+root.edge*trend,sqrt(root.edge*sigma2))       
      	for (i in N:1) x[des[i]] = rnorm(1,x[anc[i]]+el[i]*trend,sqrt(el[i]*sigma2))
      	x <- x[1:n]
      	names(x) <- phy$tip.label
      	return(x)
	}
	#----------

#------------------------------
### Simulating
if (model == "BM") {
	trans = transformed.branch.lengths(model)	
	sim = BMsim(trans, ancestral.state = p$ancestral.state, sigma2 = p$sigma2)
	}

if (model == "OU") {	
	if (p$alpha == 0) {
		trans = transformed.branch.lengths("BM")	
		sim = BMsim(trans, ancestral.state = p$ancestral.state, sigma2 = p$sigma2)
		} else sim = OUsim(y0 = p$ancestral.state, mu = p$optimal.value, alpha = p$alpha, gam = p$sigma2/(2*p$alpha))
	}

if (model=="lambda") {
	trans = transformed.branch.lengths(model, value = p$lambda)	
	sim = BMsim(trans, ancestral.state = p$ancestral.state, sigma2 = p$sigma2)
	}

if (model=="kappa") {
	trans = transformed.branch.lengths(model, value = p$kappa)	
	sim = BMsim(trans, ancestral.state = p$ancestral.state, sigma2 = p$sigma2)
	}

if (model=="delta") {
	trans = transformed.branch.lengths(model, value = p$delta)	
	sim = BMsim(trans, ancestral.state = p$ancestral.state, sigma2 = p$sigma2)
	}

if (model=="EB") {	
	trans = transformed.branch.lengths(model, value = p$rate)	
	sim = BMsim(trans, ancestral.state = p$ancestral.state, sigma2 = p$sigma2)
	}

if (model=="trend") {
	trans = transformed.branch.lengths(model)	
	sim = BMtrend.sim(trans, ancestral.state = p$ancestral.state, sigma2 = p$sigma2, trend=p$trend)
	}

return(sim)
}

