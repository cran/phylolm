\encoding{latin1}
\name{three.point.compute}
\alias{three.point.compute}
\title{Computations with a (generalized) three-point structured tree}
\description{Computes \eqn{P'V^{-1}Q} and the log(det V) of a (generalized) three-point structured matrix.}
\usage{
three.point.compute(phy, P, Q = NULL, diagMatrix = NULL, 
            check.pruningwise = TRUE, check.names = TRUE)
}
\arguments{
  \item{phy}{a rooted phylogenetic tree of type phylo with branch lengths, to represent the 3-point structured matrix \eqn{V_0}.}
  \item{P,Q}{two matrices.}
  \item{diagMatrix}{a vector containing the diagonal elements of the diagonal matrix D if V has a generalized 3-point structure.}
  \item{check.pruningwise}{the tree is assumed to be in pruningwise order if check.pruningwise = FALSE.}
  \item{check.names}{the row names of P, Q, and the names of diagMatrix are assumed to be the same as the labels of the tips in the tree if check.names = FALSE.} 
}
\value{
	\item{vec11}{\eqn{1'V^{-1}1}.}
	\item{P1}{\eqn{P'V^{-1}1}.}
	\item{PP}{\eqn{P'V^{-1}P}.}
	\item{Q1}{\eqn{Q'V^{-1}1}.}
	\item{QQ}{\eqn{Q'V^{-1}Q}.}
	\item{PQ}{\eqn{P'V^{-1}Q}.}
	\item{logd}{log(det V).}
}
\references{
Ho, L. S. T. and Ané, C. "A linear-time algorithm for Gaussian and non-Gaussian trait evolution models". Systematic Biology in review.
}
\author{Lam Si Tung Ho}
\note{
The matrix V is assumed to be \eqn{V = DV_0D} where D is the diagonal matrix with diagonal elements in "diagMatrix" and \eqn{V_0} is the 3-point structured covariance matrix determined by phy and its branch lengths.

The number of rows of P and Q and the length of diagMatrix need to be the same as the number of tips in the tree. When Q = NULL, the function only returns \eqn{1'V^{-1}1}, \eqn{P'V^{-1}1} and \eqn{P'V^{-1}P}.
}

\seealso{
\code{\link[phylolm]{transf.branch.lengths}}.
}

\examples{
tre1 = rtree(500)
tre2 = transf.branch.lengths(phy=tre1, model="OUrandomRoot",
                             parameters = list(alpha = 0.5))
x1 = rTrait(tre1)
x2 = rTrait(tre1)
y = rTrait(tre1)
one = rep(1,50)
P = cbind(one,y)
Q = cbind(x1,x2)
three.point.compute(tre2$tree,P,Q,tre2$diagMatrix)
}

