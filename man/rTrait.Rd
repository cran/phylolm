\encoding{latin1}
\name{rTrait}
\alias{rTrait}
\title{Continuous trait simulation}
\description{Simulates a continuous trait along a tree from various phylogenetic models.}
\usage{
rTrait(phy, model = c("BM", "OU", "lambda","kappa", "delta",
              "EB", "trend"), parameters = NULL)
}
\arguments{
  \item{phy}{a phylogenetic tree of type phylo with branch lengths.}
  \item{model}{a phylogenetic model. Default is "BM", for Brownian motion. Alternatives are "OU", "lambda", "kappa", "delta", "EB" and "trend".}
  \item{parameters}{List of parameters for the model (see Note).}
}
\details{
Possible phylogenetic models are the Brownian motion model (BM), the Ornstein-Uhlenbeck model (OU), Pagel's lambda model (lambda), Pagel's kappa model (kappa), Pagel's delta model (delta), the Early Burst model (EB), and the Brownian motion model with a trend (trend).
}
\value{A numeric vector with names from the tip labels in the tree.}
\author{Lam Si Tung Ho}
\note{
The default choice for the parameters are as follows: ancestral.state=0, sigma2=1, optimal.value=0 for the OU model, alpha=0 for the selection strength in the OU model, lambda=1, kappa=1, delta=1, rate=0 for the EB model, trend=0. These default choices correspond to the BM model.
}

\seealso{
\code{\link[ape]{rTraitCont}}.
}

\examples{
tre = rtree(50)
y = rTrait(phy=tre,model="OU", parameters=list(
             optimal.value=2,sigma2=1,alpha=0.1))
}


