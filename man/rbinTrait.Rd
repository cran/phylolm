\encoding{latin1}
\name{rbinTrait}
\alias{rbinTrait}
\title{Binary trait simulation}
\description{Simulates a binary trait along a phylogeny, according to the model in Ives and Garland (2010).}
\usage{
rbinTrait(phy, beta, alpha, X = NULL, model = c("LogReg"))
}
\arguments{
  \item{phy}{a phylogenetic tree of type phylo with brach lengths.}
  \item{beta}{a vector of coefficients for the logistic regression model.}
  \item{alpha}{the phylogenetic correlation parameter.}
  \item{X}{a design matrix with one row for each tip in the phylogenetic tree}
  \item{model}{Currently, only phylogenetic logistic regression is implemented.}
}

\value{a numeric vector of 0-1 values with names from the tip labels in the tree.}

\references{
Ives, A. R. and T. Garland, Jr. 2010. "Phylogenetic logistic regression for binary dependent variables". Systematic Biology \bold{59}:9-26.
}
\author{Lam Si Tung Ho}

\note{In the absence of a design matrix X, a single intercept is used, beta should be a vector of length one and the model reduces to a 2-state Markov process on the tree with stationary mean exp(beta)/(1+exp(beta)). If a design matrix X is provided, the length of beta should be equal to the number of columns in X.}

\seealso{
\code{\link[phylolm]{rTrait}}.
}
\examples{
tre = rtree(50)
x = rTrait(phy=tre)
X = cbind(rep(1,50),x)
y = rbinTrait(phy=tre, beta=c(-1,0.5), alpha=1, X=X)
}
