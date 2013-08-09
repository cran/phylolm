phylolm <- function(formula, data=list(), phy, 
	model=c("BM","OUrandomRoot","OUfixedRoot","lambda","kappa","delta","EB","trend"),
	lower.bound=NULL, upper.bound=NULL, starting.value=NULL, ...) 
{
	#----------
	### initialize	
	if (!inherits(phy, "phylo")) stop("object \"phy\" is not of class \"phylo\".")
	model = match.arg(model)	
	if ((model=="trend")&(is.ultrametric(phy)))
		stop("the trend is unidentifiable for ultrametric trees.")	
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
		warning("no tip labels, order assumed to be the same as in the tree.\n")
		data.names = phy$tip.label 
		}
	else data.names = rownames(mf)
	order = match(data.names, phy$tip.label)
	if (sum(is.na(order))>0) {
		warning("data names do not match with the tip labels.\n")
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
	transf.branch.lengths <- function(model, value = NULL)	{
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
	### a function using three point structure to do computation
	three.point.compute <- function(trans,y,X) {		
		zero <- rep(0,n + phy$Nnode)
		vec11 <- rep(0,n + phy$Nnode) # 1'V^{-1}1
		y1 <- rep(0,n + phy$Nnode) # y'V^{-1}1
		yy <- rep(0,n + phy$Nnode) # y'V^{-1}y
		logd <- rep(0,n + phy$Nnode) # logdetV
		X1 <- matrix(0,d,n + phy$Nnode) # X'V^{-1}1
		XX <- array(0,dim=c(d,d,n + phy$Nnode)) # X'V^{-1}X
		Xy <- matrix(0,d,n + phy$Nnode) # X'V^{-1}y
		root.edge = trans[[2]]
		for (i in 1:N) {
			el = trans[[1]][i]
			if (des[i] <= n) {		
				if (el > 0) {
		   			logd[des[i]] = log(el)
		   			Xy[,des[i]] = X[des[i],]*y[des[i]]/el
		   			yy[des[i]] = y[des[i]]^2/el
		   			y1[des[i]] = y[des[i]]/el
		   			XX[,,des[i]] = X[des[i],]%*%t(X[des[i],])/el
		   			X1[,des[i]] = X[des[i],]/el
		   			vec11[des[i]] = 1/el		  
					} else zero[anc[i]] = zero[anc[i]]+1				
				} else {
				if (((el<=0) & (zero[des[i]]>0)) | (zero[des[i]]>1))
					stop("logLikelihood = -Inf.")	
				logd[des[i]] = logd[des[i]] + log(1+el*vec11[des[i]])
				Xy[,des[i]] = Xy[,des[i]] - el/(1+el*vec11[des[i]])*X1[,des[i]]*y1[des[i]]
		   		yy[des[i]] = yy[des[i]] - el/(1+el*vec11[des[i]])*y1[des[i]]^2
		   		y1[des[i]] = y1[des[i]]*(1 - 1/(1+1/el/vec11[des[i]]))
				XX[,,des[i]] = XX[,,des[i]] - el/(1+el*vec11[des[i]])*X1[,des[i]]%*%t(X1[,des[i]])
		   		X1[,des[i]] = X1[,des[i]]*(1 - 1/(1+1/el/vec11[des[i]]))		   	
				vec11[des[i]] = vec11[des[i]]*(1 - 1/(1+1/el/vec11[des[i]]))
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
		return(list(vec11=vec11[ROOT],y1=y1[ROOT],yy=yy[ROOT],X1=X1[,ROOT],XX=XX[,,ROOT],Xy=Xy[,ROOT],logd=logd[ROOT])) 
		}
	#----------
	
	#----------
	### Estimate with three point stucture
	estimate.three.point <- function(trans,y,X) {
		comp = three.point.compute(trans,y,X)
		invXX = solve(comp$XX)	
		betahat = invXX%*%comp$Xy
		sigma2hat = as.numeric((comp$yy - 2*t(betahat)%*%comp$Xy + t(betahat)%*%comp$XX%*%betahat)/n)
		if (sigma2hat<0) {
			tmp = X%*%betahat - y
			sigma2hat = three.point.compute(trans,tmp,X)[[3]]/n
			}
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
			trans = transf.branch.lengths(model, value = value)
			if (flag == 1) {
				BMest = estimate.three.point(trans, y=exp(-value*D)*y, X=exp(-value*D)*X)
				return(-BMest$llh + value*sum(D)) 
				} else {
					BMest = estimate.three.point(trans, y=y, X=X)
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
			if (matchbound == 1) warning(paste("the estimation of", parameter, "matches the upper/lower bound for this parameter.\n You may change the bounds using options \"upper.bound\" and \"lower.bound\".\n",sep=" "))
			}
		trans = transf.branch.lengths(model, value = MLEvalue)
		if (flag == 1) {
			BMest = estimate.three.point(trans, y=exp(-MLEvalue*D)*y, X=exp(-MLEvalue*D)*X)
			return(list(llh=BMest$llh-MLEvalue*sum(D), betahat=BMest$betahat, sigma2hat=BMest$sigma2hat, MLEvalue=MLEvalue, vcov = BMest$vcov))
			} else {
				BMest = estimate.three.point(trans, y=y, X=X)
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
		trans = transf.branch.lengths(model)	
		BMest = estimate.three.point(trans, y=y, X=X)
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

################################################
################################################

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

################################################
################################################

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

################################################
################################################

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

################################################
################################################

residuals.phylolm <-function(object,type=c("response"), ...)
{
	type <- match.arg(type)
    	r <- object$residuals
	r	 
}

################################################
################################################

vcov.phylolm <- function(object, ...)
{
	vcov = object$vcov
	vcov
}

################################################
################################################

logLik.phylolm <- function(object, ...)
{
	res = list(logLik = object$logLik, df = object$p)
	class(res) = "logLik.phylolm"
	res
}

################################################
################################################

print.logLik.phylolm <- function (x, ...) 
{
	cat("'log Lik.' ",x$logLik," (df=",x$df,")\n", sep = "")
}

################################################
################################################

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

################################################
################################################

plot.phylolm <-function(x, ...)
{
	plot(x$y, fitted(x), xlab = "Observed value", ylab = "Fitted value")
}

################################################
################################################

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
		stop("the trend is unidentifiable for ultrametric trees.")
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
	transf.branch.lengths <- function(model, value=NULL)
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
		trans = transf.branch.lengths(model)	
		sim = BMsim(trans, ancestral.state = p$ancestral.state, sigma2 = p$sigma2)
		}

	if (model == "OU") {	
		if (p$alpha == 0) {
			trans = transf.branch.lengths("BM")	
			sim = BMsim(trans, ancestral.state = p$ancestral.state, sigma2 = p$sigma2)
			} else sim = OUsim(y0 = p$ancestral.state, mu = p$optimal.value, alpha = p$alpha, gam = p$sigma2/(2*p$alpha))
		}

	if (model=="lambda") {
		trans = transf.branch.lengths(model, value = p$lambda)	
		sim = BMsim(trans, ancestral.state = p$ancestral.state, sigma2 = p$sigma2)
		}

	if (model=="kappa") {
		trans = transf.branch.lengths(model, value = p$kappa)	
		sim = BMsim(trans, ancestral.state = p$ancestral.state, sigma2 = p$sigma2)
		}

	if (model=="delta") {
		trans = transf.branch.lengths(model, value = p$delta)	
		sim = BMsim(trans, ancestral.state = p$ancestral.state, sigma2 = p$sigma2)
		}

	if (model=="EB") {	
		trans = transf.branch.lengths(model, value = p$rate)	
		sim = BMsim(trans, ancestral.state = p$ancestral.state, sigma2 = p$sigma2)
		}

	if (model=="trend") {
		trans = transf.branch.lengths(model)	
		sim = BMtrend.sim(trans, ancestral.state = p$ancestral.state, sigma2 = p$sigma2, trend=p$trend)
		}

	return(sim)
}

################################################
################################################

transf.branch.lengths <- function(phy, 
	model = c("BM","OUrandomRoot","OUfixedRoot","lambda","kappa","delta","EB","trend"), 
	parameters = NULL, check.pruningwise = TRUE)	
{
	if (!inherits(phy, "phylo")) stop("object \"phy\" is not of class \"phylo\".")
	model = match.arg(model)	
	if ((model=="trend")&(is.ultrametric(phy)))
		stop("the trend is unidentifiable for ultrametric trees.")	
      if (is.null(phy$edge.length)) stop("the tree has no branch lengths.")
	if (is.null(phy$tip.label)) stop("the tree has no tip labels.")	
	tol = 1e-10	
	if (check.pruningwise) phy = reorder(phy,"pruningwise")
	phy2 = phy
	n <- length(phy$tip.label)
    	N <- dim(phy$edge)[1]
    	ROOT <- n + 1L
	anc <- phy$edge[, 1]
    	des <- phy$edge[, 2]

	### Default parameters
	parameters.default = c(0,1,1,1,0)
	names(parameters.default) = c("alpha", "lambda", "kappa", "delta", "rate")

	### User defined parameters
	if (is.null(parameters)) {
		parameters = parameters.default
		} else { 
			if (class(parameters)!= "list") {
				stop("please specify parameters as a list().")
				} else {
					specified <- !c(is.null(parameters$alpha),
						is.null(parameters$lambda),
						is.null(parameters$kappa),
						is.null(parameters$delta),
						is.null(parameters$rate))
					parameters.user <- c(parameters$alpha,
						parameters$lambda,
						parameters$kappa,
						parameters$delta,
						parameters$rate)
					names(parameters.default) = c("alpha", "lambda", "kappa", "delta", "rate")
					parameters <- parameters.default
		 			parameters[specified] <- parameters.user 
					}				
			}
	p = list(alpha = parameters[1],
		lambda = parameters[2],
		kappa = parameters[3],
		delta = parameters[4],
		rate = parameters[5])

	OU = c("OUrandomRoot","OUfixedRoot")
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

	### BM model
	if (model %in% c("BM","trend")) {
		edge.length = phy$edge.length
		root.edge = 0
		diagMatrix = rep(1,n)
		}
	
	### OUrandomRoot model	
	if (model=="OUrandomRoot") {
		alpha = p$alpha
		edge.length = numeric(N)
		distFromRoot <-  exp(-2*alpha*times)
		externalEdge <- des<=n
		for (i in 1:N) {
			d1 = distFromRoot[which(names(times)==anc[i])]
			if (externalEdge[i]) {d2 = exp(-2*alpha*D[des[i]])} else d2 = distFromRoot[which(names(times)==des[i])]
			edge.length[i] = d2 - d1
			}
		root.edge = min(distFromRoot)
		diagMatrix = exp(alpha*D)
		}

	### OUfixedRoot model	
	if (model=="OUfixedRoot") {
		alpha = p$alpha
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
		diagMatrix = exp(alpha*D)
		}

	### lambda model
	if (model=="lambda") {
		lambda = p$lambda
		edge.length = numeric(N)
		distFromRoot <- pruningwise.distFromRoot(phy)
		externalEdge <- des<=n
		for (i in 1:N) {				
			if (externalEdge[i]) {edge.length[i]  = phy$edge.length[i]*lambda+(1-lambda)*distFromRoot[des[i]]}
				else edge.length[i] = phy$edge.length[i]*lambda
			}
		root.edge = 0
		diagMatrix = rep(1,n)
		}
	### kappa model
	if (model=="kappa") {
		kappa = p$kappa
		edge.length = phy$edge.length^kappa
		root.edge = 0
		diagMatrix = rep(1,n)
		}

	### delta model
	if (model=="delta") {
		delta = p$delta
		edge.length = numeric(N)
		distFromRoot <- pruningwise.distFromRoot(phy)
		depth = max(distFromRoot)
		for (i in 1:N) edge.length[i] = (distFromRoot[des[i]]^delta - distFromRoot[anc[i]]^delta)*depth^(1-delta)
		root.edge = 0
		diagMatrix = rep(1,n)
		}

	### EB model
	if (model=="EB") {
		rate = p$rate
		edge.length = numeric(N)
		if (rate==0) edge.length = phy$edge.length
			else {
				distFromRoot <- pruningwise.distFromRoot(phy)
				for (i in 1:N) edge.length[i] = (exp(rate*distFromRoot[des[i]])-exp(rate*distFromRoot[anc[i]]))/rate
				}			
		root.edge = 0
		diagMatrix = rep(1,n)
		}
		
	phy2$edge.length = edge.length
	phy2$root.edge = root.edge
	names(diagMatrix) = phy2$tip.label
	return(list(tree = phy2, diagMatrix = diagMatrix))
	}

################################################
################################################

three.point.compute <- function(phy, P, Q = NULL, diagMatrix = NULL, 
	check.pruningwise = TRUE, check.names = TRUE) 
{	
	if (!inherits(phy, "phylo")) stop("object \"phy\" is not of class \"phylo\".")
	if (!check.pruningwise)	phy = reorder(phy,"pruningwise")
	if ((!check.names)&(check.pruningwise)) stop("check.names has to be TRUE unless check.pruningwise=FALSE")
	n <- length(phy$tip.label)
    	N <- dim(phy$edge)[1]
    	ROOT <- n + 1L
	anc <- phy$edge[, 1]
    	des <- phy$edge[, 2]
	if (is.null(diagMatrix)) {
		diagMatrix = rep(1,n)
		names(diagMatrix) = phy$tip.label
		}
	flag = 0
	if (is.null(Q)) {
		flag = 1
		Q = rep(1,n)
		names(Q) = phy$tip.label
		}
	P = as.matrix(P)
	Q = as.matrix(Q)	
	if (check.names) {
		if (is.null(rownames(P))) stop("P needs to have row names.")
		order = match(rownames(P), phy$tip.label)
		if (sum(is.na(order))>0) stop("row names of P does not match the names of the tips of the tree.")
		if (is.null(rownames(Q))) stop("Q needs to have row names.")
		order = match(rownames(Q), phy$tip.label)
		if (sum(is.na(order))>0) stop("row names of Q does not match the names of the tips of the tree.")
		if (is.null(names(diagMatrix))) stop("diagMatrix needs to have names.")
		order = match(names(diagMatrix), phy$tip.label)
		if (sum(is.na(order))>0) stop("names of diagMatrix does not match the names of the tips of the tree.")
		}
	if (nrow(P)!=n) stop("the number of row of P needs to be the same as the number of the tips in the tree.")
	if (nrow(Q)!=n) stop("the number of row of Q needs to be the same as the number of the tips in the tree.")
	if (length(diagMatrix)!=n) stop("the length of diagMatrix needs to be the same as the number of the tips in the tree.")
	P = cbind(rep(1,n),P)
	Q = cbind(rep(1,n),Q)
	P = P/diagMatrix
	Q = Q/diagMatrix
	colP = ncol(P)
	colQ = ncol(Q)
	zero <- rep(0,n + phy$Nnode)
	tmp11 <- rep(0,n + phy$Nnode) 
	tmpP1 <- matrix(0,colP,n + phy$Nnode)
	tmpQ1 <- matrix(0,colQ,n + phy$Nnode)
	PP <- array(0,dim=c(colP,colP,n + phy$Nnode)) 	
	logd <- rep(0,n + phy$Nnode) 
	QQ <- array(0,dim=c(colQ,colQ,n + phy$Nnode)) 
	PQ <- array(0,dim=c(colP,colQ,n + phy$Nnode)) 
	root.edge = phy$root.edge
	for (i in 1:N) {
		el = phy$edge.length[i]
		if (des[i] <= n) {		
			if (el > 0) {
	   			logd[des[i]] = log(el)
	   			PQ[,,des[i]] = P[des[i],]%*%t(Q[des[i],])/el
	   			PP[,,des[i]] = P[des[i],]%*%t(P[des[i],])/el	   		
				tmpP1[,des[i]] = P[des[i],]/el
	   			QQ[,,des[i]] = Q[des[i],]%*%t(Q[des[i],])/el
				tmpQ1[,des[i]] = Q[des[i],]/el
				tmp11[des[i]] = 1/el		  
				} else zero[anc[i]] = zero[anc[i]]+1				
			} else {
			if (((el<=0) & (zero[des[i]]>0)) | (zero[des[i]]>1))
				stop("logLikelihood = -Inf.")	
			logd[des[i]] = logd[des[i]] + log(1+el*tmp11[des[i]])
			PQ[,,des[i]] = PQ[,,des[i]] - el/(1+el*tmp11[des[i]])*tmpP1[,des[i]]%*%t(tmpQ1[,des[i]])
	   		PP[,,des[i]] = PP[,,des[i]] - el/(1+el*tmp11[des[i]])*tmpP1[,des[i]]%*%t(tmpP1[,des[i]])
			tmpP1[,des[i]] = tmpP1[,des[i]]*(1 - 1/(1+1/el/tmp11[des[i]]))		   	
			QQ[,,des[i]] = QQ[,,des[i]] - el/(1+el*tmp11[des[i]])*tmpQ1[,des[i]]%*%t(tmpQ1[,des[i]])	
			tmpQ1[,des[i]] = tmpQ1[,des[i]]*(1 - 1/(1+1/el/tmp11[des[i]]))		   	
			tmp11[des[i]] = tmp11[des[i]]*(1 - 1/(1+1/el/tmp11[des[i]]))
	   		}
		logd[anc[i]] = logd[anc[i]] + logd[des[i]]
		PQ[,,anc[i]] = PQ[,,anc[i]] + PQ[,,des[i]]
		PP[,,anc[i]] = PP[,,anc[i]] + PP[,,des[i]]
		tmpP1[,anc[i]] = tmpP1[,anc[i]] + tmpP1[,des[i]]
		QQ[,,anc[i]] = QQ[,,anc[i]] + QQ[,,des[i]]
		tmpQ1[,anc[i]] = tmpQ1[,anc[i]] + tmpQ1[,des[i]]
		tmp11[anc[i]] = tmp11[anc[i]] + tmp11[des[i]]				
	}
	logd[ROOT] = logd[ROOT] + log(1+root.edge*tmp11[ROOT]) + 2*sum(log(diagMatrix))
	PQ[,,ROOT] = PQ[,,ROOT] - root.edge/(1+root.edge*tmp11[ROOT])*tmpP1[,ROOT]%*%t(tmpQ1[,ROOT])
	PP[,,ROOT] = PP[,,ROOT] - root.edge/(1+root.edge*tmp11[ROOT])*tmpP1[,ROOT]%*%t(tmpP1[,ROOT])	
	tmpP1[,ROOT] = tmpP1[,ROOT]*(1 - 1/(1+1/root.edge/tmp11[ROOT]))
	QQ[,,ROOT] = QQ[,,ROOT] - root.edge/(1+root.edge*tmp11[ROOT])*tmpQ1[,ROOT]%*%t(tmpQ1[,ROOT])
	tmpQ1[,ROOT] = tmpQ1[,ROOT]*(1 - 1/(1+1/root.edge/tmp11[ROOT]))
	tmp11[ROOT] = tmp11[ROOT]*(1 - 1/(1+1/root.edge/tmp11[ROOT]))
	if (flag==1) return(list(vec11=PP[1,1,ROOT],P1=PP[1,-1,ROOT],PP=PP[-1,-1,ROOT],logd=logd[ROOT]))
	return(list(vec11=PP[1,1,ROOT],P1=PP[1,-1,ROOT],PP=PP[-1,-1,ROOT],Q1=QQ[1,-1,ROOT],QQ=QQ[-1,-1,ROOT],PQ=PQ[-1,-1,ROOT],logd=logd[ROOT])) 
}

################################################
################################################

phyloglm <- function(formula, data=list(), phy, model=c("LogReg"), btol = 10, log.alpha.bound = 4) 
{
	#----------
	### initialize	
	if (!inherits(phy, "phylo")) stop("object \"phy\" is not of class \"phylo\".")	
      if (is.null(phy$edge.length)) stop("the tree has no branch lengths.")
	if (is.null(phy$tip.label)) stop("the tree has no tip labels.")
	model = match.arg(model)
	phy = reorder(phy,"pruningwise")
	n <- length(phy$tip.label)
    	N <- dim(phy$edge)[1]
    	ROOT <- n + 1L
	anc <- phy$edge[, 1]
    	des <- phy$edge[, 2]
	mf = model.frame(formula=formula,data=data)
	if (nrow(mf)!=length(phy$tip.label)) stop("number of rows of the design matrix does not match with length of the tree.")
	if (is.null(rownames(mf))) {
		warning("no tip labels, order assumed to be the same as in the tree.\n")
		data.names = phy$tip.label 
		}
	else data.names = rownames(mf)
	order = match(data.names, phy$tip.label)
	if (sum(is.na(order))>0) {
		warning("data names do not match with the tip labels.\n")
		rownames(mf) = data.names
		} else {
			tmp = mf
			rownames(mf) = phy$tip.label
			mf[order,] = tmp[1:nrow(tmp),]
		}
	X = model.matrix(attr(mf, "terms"), data=mf)
 	y = model.response(mf)	
	if ( sum(!(y %in% c(0,1))) )
         stop("The model by Ives and Garland requires a binary response (dependent variable).")
	if (var(y)==0) stop("the response (dependent variable) is always 0 or always 1.")
	d = ncol(X)
	btouch = 0	
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
	### preparing for generalized tree-structure	
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
	#----------

	#----------
	transf.branch.lengths <- function(B,lL) {
		if (d> 1) g = X%*%B else g = rep(1,n)*B
		mu = as.vector(1/(1+exp(-g)))
		p = mean(mu)		
		m = rep(NA,n)		
		for (i in 1:n) 
			if (mu[i] < p) m[i] = mu[i]*sqrt(1-p)/sqrt(p) else m[i] = (1-mu[i])*sqrt(p)/sqrt(1-p)			
		alpha = 1/exp(lL)
		edge.length = numeric(N)
		distFromRoot <-  exp(-2*alpha*times)
		externalEdge <- des<=n
		for (i in 1:N) {
			d1 = distFromRoot[which(names(times)==anc[i])]
			if (externalEdge[i]) {
				if (mu[des[i]] < p) m[des[i]] = mu[des[i]]*sqrt(1-p)/sqrt(p) else m[des[i]] = (1-mu[des[i]])*sqrt(p)/sqrt(1-p)
				d2 = exp(-2*alpha*D[des[i]])*mu[des[i]]*(1-mu[des[i]])/(m[des[i]])^2				
				} else d2 = distFromRoot[which(names(times)==des[i])]
			edge.length[i] = d2 - d1			
			if (is.nan(edge.length[i])) {
				stop("edge.length[i] is NaN. Please reduce btol and/or log.alpha.bound.")
				}
			}
		root.edge = min(distFromRoot)	
		diag = m*exp(alpha*D)	
		return(list(edge.length,root.edge,diag))	
		}
	#----------

	#----------	
	### a function using three point structure to do computation
	three.point.compute <- function(trans,y,X) {	
		zero <- rep(0,n + phy$Nnode)
		y = as.vector(y)
		X = as.matrix(X)
		vec11 <- rep(0,n + phy$Nnode) # 1'V^{-1}1
		y1 <- rep(0,n + phy$Nnode) # y'V^{-1}1
		yy <- rep(0,n + phy$Nnode) # y'V^{-1}y
		logd <- rep(0,n + phy$Nnode) # logdetV
		X1 <- matrix(0,d,n + phy$Nnode) # X'V^{-1}1
		XX <- array(0,dim=c(d,d,n + phy$Nnode)) # X'V^{-1}X
		Xy <- matrix(0,d,n + phy$Nnode) # X'V^{-1}y
		root.edge = trans[[2]]
		for (i in 1:N) {
			el = trans[[1]][i]
			if (des[i] <= n) {		
				if (el > 0) {
			   		logd[des[i]] = log(el)
			   		Xy[,des[i]] = X[des[i],]*y[des[i]]/el
			   		yy[des[i]] = y[des[i]]^2/el
			   		y1[des[i]] = y[des[i]]/el
			   		XX[,,des[i]] = X[des[i],]%*%t(X[des[i],])/el
			   		X1[,des[i]] = X[des[i],]/el
			   		vec11[des[i]] = 1/el		  
					} else zero[anc[i]] = zero[anc[i]]+1				
				} else {
				if (((el<=0) & (zero[des[i]]>0)) | (zero[des[i]]>1))
					stop("logLikelihood = -Inf.")	
				logd[des[i]] = logd[des[i]] + log(1+el*vec11[des[i]])
				Xy[,des[i]] = Xy[,des[i]] - el/(1+el*vec11[des[i]])*X1[,des[i]]*y1[des[i]]
			   	yy[des[i]] = yy[des[i]] - el/(1+el*vec11[des[i]])*y1[des[i]]^2
			   	y1[des[i]] = y1[des[i]]*(1 - 1/(1+1/el/vec11[des[i]]))
				XX[,,des[i]] = XX[,,des[i]] - el/(1+el*vec11[des[i]])*X1[,des[i]]%*%t(X1[,des[i]])
			   	X1[,des[i]] = X1[,des[i]]*(1 - 1/(1+1/el/vec11[des[i]]))		   	
				vec11[des[i]] = vec11[des[i]]*(1 - 1/(1+1/el/vec11[des[i]]))
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
		return(list(vec11=vec11[ROOT],y1=y1[ROOT],yy=yy[ROOT],X1=X1[,ROOT],XX=XX[,,ROOT],Xy=Xy[,ROOT],logd=logd[ROOT])) 
		}
	#----------

	#----------
	plogregfunct <- function(startB,startlL) {
		convergeflag = 1
		clL = startlL
		diflL = 100
		difB = 100
		
		counter = 0
		cB = startB

		ttozero = 10^6
		while (((diflL>10^-6)|(difB>10^-6)|(ttozero>10^-1))&(counter<10)) {
			counter = counter+1 

			oldlL = clL
			oldB = cB
			
			olddiflL = diflL
			olddifB = difB

			opt <- optim(par = cB, fn = function(par){plogregBfunct(par,clL)}, method = "L-BFGS-B", control = list(factr=1e12))
			cB = as.vector(opt$par)
			ttozero = as.numeric(opt$value)
			difB = t(cB-oldB)%*%(cB-oldB)
			
			if (ttozero > 10^-2) {				
				Btemp = btol/5*rnorm(length(startB),0,1) + startB
				opt <- optim(par = Btemp, fn = function(par){plogregBfunct(par,clL)}, method = "L-BFGS-B", control = list(factr=1e12))
				Btemp = as.vector(opt$par)
				newttozero = as.numeric(opt$value)
				if (newttozero < ttozero) {
					cB = Btemp
					ttozero = newttozero 
					}
				difB = t(cB-oldB)%*%(cB-oldB)				
				} 
			opt <- optim(par = clL, fn = function(par){plogreglLfunct(cB,par)}, method = "L-BFGS-B")
			clL = as.numeric(opt$par)
			diflL = (clL-oldlL)^2
			}
			
		if (counter>=10) {
			clL = (clL+oldlL)/2
			cB = (cB+oldB)/2
			while (((diflL>10^-6)|(difB>10^-6)|(ttozero>10^-1))&(counter<20)) {
				counter = counter + 1
				
				clL = (clL+oldlL)/2
				oldB = cB
				
				olddiflL = diflL
				olddifB = difB
			
				opt <- optim(par = cB, fn = function(par){plogregBfunct(par,clL)}, method = "L-BFGS-B", control = list(factr=1e12))
				cB = as.vector(opt$par)
				ttozero = as.numeric(opt$value)
				difB = t(cB-oldB)%*%(cB-oldB)
				if ((ttozero > 10^-2)&(abs(clL-log(Tmax))<log.alpha.bound)) {				
					Btemp = btol/5*rnorm(length(startB),0,1) + startB
					opt <- optim(par = Btemp, fn = function(par){plogregBfunct(par,clL)}, method = "L-BFGS-B", control = list(factr=1e12))
					Btemp = as.vector(opt$par)
					newttozero = as.numeric(opt$value)
					if (newttozero < ttozero) {
						cB = Btemp
						ttozero = newttozero 
						}
					difB = t(cB-oldB)%*%(cB-oldB)				
					} 				
				cB = (cB+oldB)/2
				oldlL = clL

				opt <- optim(par = clL, fn = function(par){plogreglLfunct(cB,par)}, method = "L-BFGS-B")
				clL = as.numeric(opt$par)
				diflL = (clL-oldlL)^2
				}	
			if (counter >= 19) 
				if ((max(abs(c(oldlL - clL,oldB - cB)))>0.1)|(ttozero > 0.5)) convergeflag = 0					
			}	
		return(list(B=cB,lL=clL,convergeflag = convergeflag))	
		}
	#----------	

	#----------	
	plogreg1Dfunct <- function(startB,startlL) {
		convergeflag = 1
		clL = startlL
		diflL = 100
		difB = 100
		
		counter = 0
		cB = startB

		ttozero = 10^6
		while (((diflL>10^-3)|(difB>10^-3)|(ttozero>10^-1))&(counter<20)) {
			counter = counter+1 

			oldlL = clL
			oldB = cB
			
			olddiflL = diflL
			olddifB = difB

			opt <- optim(par = cB, fn = function(par){plogregB1Dfunct(par,clL)}, method = "L-BFGS-B")
			cB = as.numeric(opt$par)
			ttozero = as.numeric(opt$value)
			difB = (cB-oldB)^2
						
			if ((ttozero > 10^-2)&(abs(clL-log(Tmax))<log.alpha.bound)) {				
					Btemp = btol/5*rnorm(length(startB),0,1) + startB
					opt <- optim(par = Btemp, fn = function(par){plogregBfunct(par,clL)}, method = "L-BFGS-B")
					Btemp = as.vector(opt$par)
					newttozero = as.numeric(opt$value)
					if (newttozero < ttozero) {
						cB = Btemp
						ttozero = newttozero 
						}
					difB = t(cB-oldB)%*%(cB-oldB)				
					} 			

			opt <- optim(par = clL, fn = function(par){plogreglLfunct(cB,par)}, method = "L-BFGS-B")
			clL = as.numeric(opt$par)
			diflL = (clL-oldlL)^2
			}

		if (counter>=20) {
			clL = (clL+oldlL)/2
			cB = (cB+oldB)/2
			while (((diflL>10^-3)|(difB>10^-3)|(ttozero>10^-1))&(counter<40)) {
				counter = counter + 1
				
				clL = (clL+oldlL)/2
				oldB = cB
				
				olddiflL = diflL
				olddifB = difB
			
				opt <- optim(par = cB, fn = function(par){plogregB1Dfunct(par,clL)}, method = "L-BFGS-B")
				cB = as.numeric(opt$par)
				ttozero = as.numeric(opt$value)
				difB = (cB-oldB)^2
				if ((ttozero > 10^-2)&(abs(clL-log(Tmax))<log.alpha.bound)) {				
					Btemp = btol/5*rnorm(length(startB),0,1) + startB
					opt <- optim(par = Btemp, fn = function(par){plogregBfunct(par,clL)}, method = "L-BFGS-B")
					Btemp = as.vector(opt$par)
					newttozero = as.numeric(opt$value)
					if (newttozero < ttozero) {
						cB = Btemp
						ttozero = newttozero 
						}
					difB = t(cB-oldB)%*%(cB-oldB)				
					} 					
				cB = (cB+oldB)/2
				oldlL = clL

				opt <- optim(par = clL, fn = function(par){plogreglLfunct(cB,par)},method = "L-BFGS-B")
				clL = as.numeric(opt$par)
				diflL = (clL-oldlL)^2
				}	
			if (counter >= 39) 
				if (max(abs(c(oldlL - clL,oldB - cB)))>0.5) convergeflag = 0					
			}	
		return(list(B=cB,lL=clL,convergeflag = convergeflag))	
		}
	#----------	
	
	#----------
	plogregBfunct <- function(B,lL) {	
		g = X%*%B
		if (sum(abs(g) >= btol)>0) {
			btouch <<- 1
			return(1e6)
			}
		mu = as.vector(1/(1+exp(-g)))
				
		temp = transf.branch.lengths(B,lL)
		trans = temp[1:2]
		diag = temp[[3]]
		comp = three.point.compute(trans,(y-mu)/diag,mu*(1-mu)*X/diag)
		
		logdetC = comp$logd + 2*sum(log(diag)) - sum(log(mu*(1-mu)))	
		if (logdetC < -100*log(10)) return(1e6)
		
		Z = comp$Xy
		
		# Firth correction
		Dx = 0.1
		k = length(B)
		infoM = comp$XX
		invInfoM = solve(infoM)
		FirthC = rep(NA,k)
		for (i in 1:k) {
			# increase
			dB = B
			dB[i] = dB[i]+Dx
			
			g = X%*%dB
			if (sum(abs(g) >= btol)>0) return(1e6)
			mu = as.vector(1/(1+exp(-g)))
			
			ttemp = transf.branch.lengths(dB,lL)
			ttrans = ttemp[1:2]
			tdiag = ttemp[[3]]
			tcomp = three.point.compute(ttrans,(y-mu)/tdiag,mu*(1-mu)*X/tdiag)
			dinfoMp = tcomp$XX
				
			# decrease
			dB = B
			dB[i] = dB[i]-Dx
		
			g = X%*%dB
			if (sum(abs(g) >= btol)>0) return(1e6)
			mu = as.vector(1/(1+exp(-g)))
		
			ttemp = transf.branch.lengths(dB,lL)
			ttrans = ttemp[1:2]
			tdiag = ttemp[[3]]
			tcomp = three.point.compute(ttrans,(y-mu)/tdiag,mu*(1-mu)*X/tdiag)
			dinfoMm = tcomp$XX

			DinfoM = (dinfoMp - dinfoMm)/Dx/2
			FirthC[i] = sum(diag(invInfoM%*%DinfoM))/2
			}
		tozero = Z + FirthC
		return(sum(tozero^2))
		}		
	#----------

	#----------
	plogregB1Dfunct <- function(B,lL) {	
		g = rep(1,n)*B
		if (sum(abs(g) >= btol)>0) {
			btouch <<- 1
			return(1e6)
			}
		mu = as.vector(1/(1+exp(-g)))

		temp = transf.branch.lengths(B,lL)
		trans = temp[1:2]
		diag = temp[[3]]
		comp = three.point.compute(trans,(y-mu)/diag,mu*(1-mu)/diag)
		Z = comp$Xy
		FirthC = (1-2*mu)/2
		tozero = Z + FirthC
		return(sum(tozero^2))
		}
	#----------
	
	#----------
	plogreglLfunct <- function(B,lL) {
		g = X%*%B
		mu = as.vector(1/(1+exp(-g)))
	
		if (abs(lL - log(Tmax)) >= log.alpha.bound) return(1e10)
		temp = transf.branch.lengths(B,lL)
		trans = temp[1:2]
		diag = temp[[3]]
		comp = three.point.compute(trans,(y-mu)/diag,mu*(1-mu)*X/diag)
		LL = (comp$logd + 2*sum(log(diag)) + comp$yy)/2
		if (!is.finite(LL)) LL = 1e10
		return(LL)
		}	
	#----------
	
	#----------
	plogregBSEfunct <- function(B,lL) {
		g = X%*%B
		mu = as.vector(1/(1+exp(-g)))
		
		temp = transf.branch.lengths(B,lL)
		trans = temp[1:2]
		diag = temp[[3]]
		comp = three.point.compute(trans,(y-mu)/diag,mu*(1-mu)*X/diag)
		infoM = comp$XX
		Q = solve(infoM)
		BSE = sqrt(diag(Q))
		covBSE = Q
		return(list(BSE = BSE, covBSE = covBSE))
		}
	#----------
	
	#------------------------------
	### Default starting values
	### B: coefficients
	### lL = log(1/alpha), alpha: phylogenetic signal	
	fit = glm(y~X-1,family=binomial)
	startB = fit$coefficients
	if (sum(abs(X%*%startB) >= btol)>0) 
		stop("Increase the bound 'btol' for the linear predictor. \n The current value excludes the estimated coefficients in the absence of phylogenetic signal.")
	startlL = log(Tmax) 

	### Fitting
	if (d > 1) {
		plogreg = plogregfunct(startB,startlL)
		lL = plogreg$lL
		B = plogreg$B
		convergeflag = plogreg$convergeflag
	
		plogregBSE = plogregBSEfunct(B,lL)
		BSE = plogregBSE$BSE
		covBSE = plogregBSE$covBSE					
		} else {
			plogreg = plogreg1Dfunct(startB,startlL)
			lL = plogreg$lL
			B = plogreg$B
			convergeflag = plogreg$convergeflag

			plogregBSE = plogregBSEfunct(B,lL)
			BSE = plogregBSE$BSE
			covBSE = plogregBSE$covBSE		
			}

	if ((lL - log(Tmax) + 0.02) > log.alpha.bound) {
		warn = paste("the estimate of 'alpha' (",1/exp(lL),") reaches the lower bound (",1/Tmax/exp(log.alpha.bound),").\n This may reflect a flat likelihood at low alpha values near 0, meaning that the phylogenetic correlation is estimated to be maximal under the model in Ives and Garland (2010).", sep="")
		warning(warn)	
		}
	if ((lL - log(Tmax) - 0.02) < - log.alpha.bound) {
		warn = paste("the estimate of 'alpha' (",1/exp(lL),") reaches the upper bound (",exp(log.alpha.bound)/Tmax,").\n This may simply reflect a flat likelihood at large alpha values, meaning that the phylogenetic correlation is estimated to be negligible.",sep="")
		warning(warn)	
		}
	if (btouch == 1) 
		warning("the boundary of the linear predictor has been reached during the optimization procedure. \n You can increase this bound by increasing 'btol'.")
	if (convergeflag == 0) warning("phyloglm failed to converge.\n")
	results = list(coefficients = B, alpha = 1/exp(lL), sd = BSE, vcov = covBSE, convergeflag = convergeflag)
	names(results$coefficients) = colnames(X)
	colnames(results$vcov) = colnames(X)
	rownames(results$vcov) = colnames(X)
	results$linear.predictors = as.vector(X %*% results$coefficients)
	names(results$linear.predictors) = names(y)
	if (max(abs(results$linear.predictors)) + 0.01 > btol) warning("the linear predictor reaches its bound.")
	results$fitted.values = as.vector(1/(1+exp(-results$linear.predictors)))
	names(results$fitted.values ) = names(y)
	results$residuals = y - results$fitted.values
	results$mean.tip.height = Tmax
	results$y = y
	results$n = n
	results$d = d
	results$formula = formula
	results$call = match.call()
	results$model = model
	class(results) = "phyloglm"
	results
}

################################################
################################################

print.phyloglm <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
	cat("Call:\n")
	print(x$call)
	cat("\n")
	cat("\nParameter estimate(s) from GEE approximation:\n")
	cat("alpha:",x$alpha)
	cat("\n")
	
	cat("\nCoefficients:\n")
	print(x$coefficients)
}

################################################
################################################

summary.phyloglm <- function(object, ...) 
{
	se <- object$sd
	TAB <- cbind(Estimate = coef(object),
		StdErr = se)
	res <- list(call=object$call,
		coefficients=TAB,
		residuals = object$residuals,
		alpha=object$alpha,		
		model=object$model)
	class(res) = "summary.phyloglm"
	res
}

################################################
################################################

print.summary.phyloglm <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
	cat("\nCall:\n")
	print(x$call)
	cat("\n")
	

	cat("\nParameter estimate(s) from GEE approximation:\n")
	cat("alpha:",x$alpha)
	cat("\n")

	cat("\nCoefficients:\n")
	printCoefmat(x$coefficients, P.values=FALSE, has.Pvalue=FALSE)
}

################################################
################################################

residuals.phyloglm <-function(object,type=c("response"), ...)
{
	type <- match.arg(type)
    	r <- object$residuals
	r	 
}

################################################
################################################

plot.phyloglm <-function(x, ...)
{
	plot(fitted(x), x$residuals, xlab = "Fitted value", ylab = "Residuals")
	abline(h=0, lty = 2)
}

################################################
################################################

vcov.phyloglm <- function(object, ...)
{
	vcov = object$vcov
	vcov
}

################################################
################################################

rbinTrait <-function(phy, beta, alpha, X = NULL, model = c("LogReg")) {
	if (is.null(phy$edge.length)) stop("the tree has no branch lengths.")
	if (is.null(phy$tip.label)) stop("the tree has no tip labels.")		
	model = match.arg(model)
	phy = reorder(phy,"pruningwise")
	n <- length(phy$tip.label)
   	N <- dim(phy$edge)[1]
    	ROOT <- n + 1L
	anc <- phy$edge[, 1]
    	des <- phy$edge[, 2]
	B = beta
	d = length(B)	

	if ((d==1)&(!is.null(X))) stop("The design matrix is not needed when the coefficient is a scalar.")

	if (d>1) {
		if (is.null(X)) stop("there is no independent variables.")
		X = as.matrix(X)
		if (length(B)!=ncol(X)) stop("number of columns of the design matrix does not match the length of the vector of coefficients.")
		if (nrow(X)!=n) stop("number of rows of the design matrix does not match the number of tips in the tree.")
		if (is.null(rownames(X))) {
			warning("independent variables have no tip labels, order assumed to be the same as in the tree.\n")
			data.names = phy$tip.label 
			} else data.names = rownames(X)
		order = match(data.names, phy$tip.label)
		if (sum(is.na(order))>0) stop("data names do not match with the tip labels.\n")
		g = X%*%B 
		mu = as.vector(exp(g)/(1+exp(g)))
		p = mean(mu)
		} else {
			g = as.numeric(B)
			p = exp(g)/(1+exp(g))		
			}
	q = 1-p
	
	y <- numeric(n + phy$Nnode)	
	el = phy$edge.length	
	y[ROOT] = as.numeric(runif(1)<p)
      for (i in N:1)
		if (y[anc[i]] == 0) y[des[i]] = as.numeric(runif(1)<(p-p*exp(-el[i]*alpha)))
			else y[des[i]] = as.numeric(runif(1)<(p+q*exp(-el[i]*alpha)))	
      y <- y[1:n]

	if (d>1) {
		for (i in 1:n) {
			et = min(mu[i]/p,(1-mu[i])/(1-p))
			b = as.numeric(mu[i]>p)
			if (y[i] == 0) y[i] = as.numeric(runif(1)<(b-b*et)) 
				else y[i] = as.numeric(runif(1)<(b+(1-b)*et)) 
			}
		}
      names(y) <- phy$tip.label
      return(y)
}


