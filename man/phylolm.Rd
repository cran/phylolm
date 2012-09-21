\encoding{latin1}
\name{phylolm}
\alias{phylolm}

\title{Phylogenetic Linear Model}
\description{Fits a phylogenetic linear regression model. The likelihood is calculated with an algorithm that is linear in the number of tips in the tree.}
\usage{
phylolm(formula, data = list(), phy, model = c("BM", "OUrandomRoot",
        "OUfixedRoot", "lambda", "kappa", "delta", "EB", "trend"),
        lower.bound = NULL, upper.bound = NULL, starting.value = NULL, ...)
}

\arguments{
  \item{formula}{a model formula.}
  \item{data}{a data frame containing variables in the model. If not found in data, the variables are taken from current environment.}
  \item{phy}{a phylogenetic tree of type phylo.}
  \item{model}{a model for the covariance (see Details).}
  \item{lower.bound}{optional lower bound for the optimization of the phylogenetic model parameter.}
  \item{upper.bound}{optional upper bound for the optimization of the phylogenetic model parameter.}
  \item{starting.value}{optional starting value for the optimization of the phylogenetic model parameter.}
  \item{\dots}{further arguments to be passed to the function \code{optim}.}
}
\details{This function uses an algorithm that is linear in the number of tips in the tree to calculate the likelihood. Possible phylogenetic models for the error term are the Brownian motion model (BM), the Ornstein-Ühlenbeck model with an ancestral state to be estimated at the root (OUfixedRoot), the Ornstein-Ühlenbeck model with the ancestral state at the root having the stationary distribution (OUrandomRoot), Pagel's lambda model (lambda), Pagel's kappa model (kappa), Pagel's delta model (delta), the Early Burst model (EB), and the Brownian motion model with a trend (trend).}
\value{
   \item{coefficients}{the named vector of coefficients.}
   \item{sigma2}{the maximum likelihood estimate of the variance rate sigma square.}
   \item{optpar}{the optimized value of the phylogenetic correlation parameter (alpha, lambda, kappa, delta or rate).}
   \item{logLik}{the maximum of the log likelihood.}
   \item{p}{the number of all parameters of the model.}
   \item{aic}{AIC value of the model.}
   \item{vcov}{covariance matrix for the coefficients, given the phylogenetic correlation parameter (if any).}
   \item{fitted.values}{fitted values.}
   \item{residuals}{raw residuals.}
   \item{y}{response.}
   \item{n}{number of observations (tips in the tree).}
   \item{d}{number of dependent variables.}
   \item{formula}{the model formula.}
   \item{call}{the original call to the function.}
   \item{model}{the phylogenetic model for the covariance.}
}
\references{
Ané, C. and Ho, L. S. T. \emph{A linear-time algorithm for calculating the likelihood under Brownian-structured models for quantitative trait evolution}. In review.

Butler, M. A. and King, A. A. 2004. \emph{Phylogenetic comparative analysis: A modeling approach for adaptive evolution}. The American Naturalist \bold{164}:683-695.

Hansen, T. F. 1997. \emph{Stabilizing selection and the comparative analysis of adaptation}. Evolution \bold{51}:1341-1351.

Harmon, L. J. et al. 2010. \emph{Early bursts of body size and shape evolution are rare in comparative data}. Evolution \bold{64}:2385-2396.

Ho, L. S. T. and Ané, C. \emph{Asymptotic theory with hierarchical autocorrelation: Ornstein-Ühlenbeck tree models}. In review.

Pagel, M. 1997. \emph{Inferring evolutionary processes from phylogenies}. Zoologica Scripta \bold{26}:331-348.

Pagel, M. 1999. \emph{Inferring the historical patterns of biological evolution}. Nature \bold{401}:877-884.

}
\author{Lam Si Tung Ho, Cécile Ané}
\note{
The tip labels in the tree are matched to the data names (row names in the data frame). If no data frame is provided through the argument data, the taxa are assumed to come in the same order for all variables in the model. In that case, taxon labels in the tree are matched to taxon labels in the response variable based on the row names of the response vector.

For the delta model, the tree is rescaled back to its original height after each node's distance from the root is raised to the power delta. This is to provide a stable estimate of the variance parameter sigma2. For non-ultrametric trees, the tree is rescaled to keep the longest distance from the root  to its original value.

The trend model can only be used with non-ultrametric trees. For this model, one predictor variable is added to the model whose values are the distances from the root to every tip of the tree. The estimate of the coefficent for this variable forms the trend value.
}
\seealso{
\code{\link[ape]{corBrownian}}, \code{\link[ape]{corMartins}}, \code{\link[ape]{corPagel}}, \code{\link[geiger]{fitContinuous}}, \code{\link[caper]{pgls}}.
}

\examples{
p = rcoal(50)
taxa = sort(p$tip.label)
b0=0; b1=1;
x = rTrait(phy=p,model="BM",
             parameters=list(ancestral.state=0,sigma2=10))
y  = b0 + b1*x + 
     rTrait(phy=p,model="lambda",parameters=list(
              ancestral.state=0,sigma2=1,lambda=0.5))
dat = data.frame(trait=y[taxa],pred=x[taxa])
fit = phylolm(trait~pred,data=dat,phy=p,model="lambda")
summary(fit)
}


