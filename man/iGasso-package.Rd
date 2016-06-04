\name{iGasso-package}
\alias{iGasso-package}
\alias{iGasso}
\docType{package}
\title{
Statistical Tests and utilities for Genetic Association
}
\description{
\code{iGasso} is a collection of statistical tests developed by our group for genetic association studies. So far it contains functions for rare variants association, for association with multiple phenotypes, for linear mixed model analysis, and for model-free association analysis. There is also a function for genome plot. It will keep growing as more tests are developed. Use \code{?iGasso} to see an introduction.
}
\details{
\tabular{ll}{
Package: \tab iGasso\cr
Type: \tab Package\cr
Version: \tab 1.4\cr
Date: \tab 2016-06-3\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}
}
\author{
Kai Wang \code{<kai-wang@uiowa.edu>}
}
\references{
Anscombe F.J. (1948) The transformation of Poisson, binomial and negative-binomial data. \emph{Biometrika} \bold{35(3/4)}, 246--254.

Chanter, D. O. (1975). Modifications of the angular transformation. \emph{Journal of the Royal Statistical Society. Series B (Applied Statistics)}, \bold{24 (3)}, 354--359.

Freeman, M. F., Tukey, J. W. (1950) Transformations related to the angular and the square root. \emph{The Annals of Mathematical Statistics} \bold{21(4)}, 607--611.

Wang, K. (2012) An application of the proportional odds model to genetic association studies. Submitted.

Wang K. (2012) Statistical tests of genetic association for case-control study designs. \emph{Biostatistics}. Accepted. PMID: 22389176

Wang, K., Fingert, J. (2012) Statistical tests for detecting rare variants using variance-stabilizing transformations. \emph{Annals of Human Genetics}. Accepted. 

Zar, J. H. (1999) \emph{Biostatistical Analysis, 4th ed.}, New Jersey:Prentice-Hall, Inc.
}

%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
%%\keyword{ package }
%\seealso{
%%~~ Optional links to other man pages, e.g. ~~
%%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
y = rnorm(100)
chr = c(rep(1, 20), rep(3, 20), rep(10, 20), rep(19, 30), rep("X", 10))
pos = c(1:20, 1:20, 1:20, 1:30, 1:10)
mydata = data.frame(y=y, chr=chr, pos=pos)
genome.plot(mydata, sig.line=c(1, -1), ylab="T Statistic")

G = rbind(c(14, 999), c(3, 1081))
VSTF.test(G)

G = rbind(c(161, 474, 489), c(231, 444, 380))
MFree.test(G)

G = matrix(sample(c(0,1,2), 200, replace=TRUE), ncol=10)
y = rnorm(10)
X = matrix(rnorm(10), ncol=1)
}
