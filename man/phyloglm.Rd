\encoding{latin1}
\name{phyloglm}
\alias{phyloglm}

\title{Phylogenetic Generalized Linear Model}
\description{Fits the phylogenetic logistic regression described in
  Ives and Garland (2010). The computation uses an algorithm that is
  linear in the number of tips in the tree.}
\usage{
phyloglm(formula, data, phy, method = c("MPLE","IG10"), btol = 10,
         log.alpha.bound = 4, start.beta=NULL, start.alpha=NULL)
}

\arguments{
  \item{formula}{a model formula.}
  \item{data}{a data frame containing variables in the model. If not
  found in \code{data}, the variables are taken from the current environment.}
  \item{phy}{a phylogenetic tree of type phylo with branch lengths.}
  \item{method}{The "IG10" method optimizes a GEE approximation to the
  penalized likelihood. "MPLE" maximizes the penalized likelihood. In
  both cases, the penalty is Firth's correction.}
  \item{btol}{bound on the linear predictor to bound the searching space.}
  \item{log.alpha.bound}{bound for the log of the parameter alpha}
  \item{start.beta}{starting values for beta coefficients}
  \item{start.alpha}{starting values for alpha (phylogenetic correlation)}
}
\details{
This function uses an algorithm that is linear in the number of tips in the tree.
}
\value{
  \item{coefficients}{the named vector of coefficients.}
  \item{alpha}{the phylogenetic correlation parameter.}
  \item{sd}{standard deviation for the regression coefficients, given
  the phylogenetic correlation parameter \eqn{\alpha}{alpha}.}
  \item{vcov}{covariance matrix for the regression coefficients, given
    the phylogenetic correlation parameter.}
  \item{logLik}{log likelihood}
  \item{aic}{AIC}
  \item{penlogLik}{penalized log likelihood, using Firth's penalty for coefficients}
  \item{y}{response}
  \item{n}{number of observations (tips in the tree)}
  \item{d}{number of dependent variables}
  \item{formula}{the model formula}
  \item{call}{the original call to the function}
  \item{method}{the estimation method}
  \item{convergence}{An integer code with '0' for successful optimization.
  With MPLE, this is the convergence code from the \code{optim} routine.}
}
\references{
Ho, L. S. T. and Ané, C. "A linear-time algorithm for Gaussian and
  non-Gaussian trait evolution models".  Systematic Biology \bold{63}(3):397-408.

Ives, A. R. and T. Garland, Jr. 2010. "Phylogenetic logistic regression
for binary dependent variables".  Systematic Biology \bold{59}:9-26.
}
\author{Lam Si Tung Ho, Robert Lachlan and Cécile Ané}
\note{
The tip labels in the tree are matched to the data names (row names in
the data frame). If no data frame is provided through the argument \code{data},
taxon labels in the tree are matched to taxon labels in the response
variable based on the row names of the response vector, and the taxa are
assumed to come in the same order for all variables in the model.

The default choice \code{btol=10} constrains the fitted values, i.e. the
probability of "1" predicted by the model, to lie within
\eqn{1/(1+e^{ 10})=0.000045}{1/(1+exp(+10)) = 0.000045} and
\eqn{1/(1+e^{-10})=0.999955}{1/(1+exp(-10)) = 0.999955}.

The log of \eqn{\alpha}{alpha} is bounded in the interval
\eqn{-\log(T) \pm \mathrm{log.alpha.bound}}{(-log(T) - log.alpha.bound,
  -log(T) + log.alpha.bound)} 
where \eqn{T} is the mean of the distances from the root to tips. In
other words, \eqn{\alpha T}{alpha*T} is constrained to lie within
\eqn{\exp(\pm\mathrm{log.alpha.bound})}{exp(+/- log.alpha.bound)}.
}

\examples{
set.seed(123456)
tre = rtree(50)
x = rTrait(n=1,phy=tre)
X = cbind(rep(1,50),x)
y = rbinTrait(n=1,phy=tre, beta=c(-1,0.5), alpha=1 ,X=X)
dat = data.frame(trait01 = y, predictor = x)
fit = phyloglm(trait01~predictor,phy=tre,data=dat)
summary(fit)
coef(fit)
vcov(fit)
}
