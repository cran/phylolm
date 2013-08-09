\encoding{latin1}
\name{transf.branch.lengths}
\alias{transf.branch.lengths}
\title{Creates a tree with branch lengths to represent the 3-point structure of a covariance matrix}
\description{Creates a phylogenetic tree with branch lengths and a diagonal matrix to represent a (generalized) 3-point structured covariance matrix from a trait evolution model on a known tree.}
\usage{
transf.branch.lengths(phy, model = c("BM", "OUrandomRoot",
       "OUfixedRoot", "lambda", "kappa", "delta", "EB", "trend"),
       parameters = NULL, check.pruningwise = TRUE)
}
\arguments{
  \item{phy}{a phylogenetic tree of type phylo with branch lengths.}
  \item{model}{a phylogenetic model. Default is "BM", for Brownian motion. Alternatives are "OU", "lambda", "kappa", "delta", "EB" and "trend".}
  \item{parameters}{List of parameters for the model (see Note).}
  \item{check.pruningwise}{the tree is assumed to be in pruningwise order if check.pruningwise = FALSE.}
}
\details{
Possible phylogenetic models are the Brownian motion model (BM), the Ornstein-Uhlenbeck model (OU), Pagel's lambda model (lambda), Pagel's kappa model (kappa), Pagel's delta model (delta), the Early Burst model (EB), and the Brownian motion model with a trend (trend).}
\value{
	\item{tree}{a rooted tree with a root edge and transformed branch lengths.}
	\item{diagMatrix}{a vector containing the diagonal elements of the diagonal matrix D if the covariance matrix has a generalized 3-point structure.}
}
\references{
Ho, L. S. T. and Ané, C. \emph{A linear-time algorithm for Gaussian and non-Gaussian trait evolution models}. Systematic Biology in review.
}
\author{Lam Si Tung Ho}
\note{
The default choice for the parameters are as follows: alpha=0 for the selection strength in the OU model, lambda=1, kappa=1, delta=1, rate=0 for the EB model. These default choices correspond to the BM model.

Edges in the output tree are in pruningwise order.

If model = "BM" or model ="trend", the output tree is the same as the input tree except that the output tree is in pruningwise order.
}

\seealso{
\code{\link[phylolm]{three.point.compute}}.
}

\examples{
set.seed(123456)
tre1 = rcoal(10)
tre2 = transf.branch.lengths(phy=tre1, model="OUrandomRoot",
                             parameters = list(alpha=1))
par(mfrow = c(2,1))
plot(tre1)
plot(tre2$tree,root.edge=TRUE)
}
