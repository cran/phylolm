\encoding{latin1}
\name{phyloglm}
\alias{phyloglm}

\title{Phylogenetic Generalized Linear Model}
\description{Fits a phylogenetic generalized linear model. The computation uses an algorithm that is linear in the number of tips in the tree.}
\usage{
phyloglm(formula, data = list(), phy, model = c("LogReg"), 
                  btol = 10, log.alpha.bound = 4)
}

\arguments{
  \item{formula}{a model formula.}
  \item{data}{a data frame containing variables in the model. If not found in data, the variables are taken from the current environment.}
  \item{phy}{a phylogenetic tree of type phylo with branch lengths.}
  \item{model}{Currently, only phylogenetic logistic regression is implemented.}
  \item{btol}{bound on the linear predictor to bound the searching space.}
  \item{log.alpha.bound}{bound for the log of the parameter alpha.}
}
\details{
This function uses an algorithm that is linear in the number of tips in the tree.
}
\value{
	\item{coefficients}{the named vector of coefficients.}
	\item{alpha}{the phylogenetic correlation parameter.}
	\item{sd}{standard deviation for the regression coefficients, given the phylogenetic correlation parameter.}
	\item{vcov}{covariance matrix for the regression coefficients, given the phylogenetic correlation parameter.}
	\item{y}{response.}
	\item{n}{number of observations (tips in the tree).}
      \item{d}{number of dependent variables.}
  	\item{formula}{the model formula.}
   	\item{call}{the original call to the function.}
   	\item{model}{the phylogenetic model for the covariance.}
}
\references{
Ho, L. S. T. and Ané, C. "A linear-time algorithm for Gaussian and non-Gaussian trait evolution models". Systematic Biology in review.

Ives, A. R. and T. Garland, Jr. 2010. "Phylogenetic logistic regression for binary dependent variables". Systematic Biology \bold{59}:9-26.
}
\author{Lam Si Tung Ho}
\note{
The tip labels in the tree are matched to the data names (row names in the data frame). If no data frame is provided through the argument data, taxon labels in the tree are matched to taxon labels in the response variable based on the row names of the response vector, and the taxa are assumed to come in the same order for all variables in the model.

The default choice btol=10 constrains the fitted values, i.e. the probability of "1" predicted by the model, to lie within 1/(1+exp(10))=0.000045 and 1/(1+exp(-10)) = 0.999955.

The log of the parameter alpha is bounded in the interval (-log(d) - log.alpha.bound, -log(d) + log.alpha.bound) where d is the mean of the distances from the root to tips.
}

\examples{
set.seed(123456)
tre = rtree(50)
x = rTrait(phy=tre)
X = cbind(rep(1,50),x)
y = rbinTrait(phy=tre, beta=c(-1,0.5), alpha=1 ,X=X)
dat = data.frame(trait01 = y, predictor = x)
fit = phyloglm(trait01~predictor,phy=tre,data=dat)
summary(fit)
coef(fit)
vcov(fit)
}
