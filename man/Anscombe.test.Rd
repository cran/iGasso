\name{Anscombe.test}
\alias{Anscombe.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genetic Association Test Based on the Anscombe Transformation.
}
\description{
\code{Anscombe.test} performs a test on association between a variant and case-control status using the Anscombe transformation.
}
\usage{
Anscombe.test(G)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G }{
a \code{2 x 2} matrix. The first row is for cases and the second one for controls. In each row, the first element is the number of non-carriers and the second one is the number of carriers with at least 1 copy of the variant. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with class "\code{htest}" containing the following components:
\item{statistic}{the value of the test statistic.}
\item{p.value}{the p-value for the test computed from a chi-square distribution with 1 df.}
\item{method}{a character string indicting the test performed.}
\item{data.name}{a character string giving the name of the data.}
}
\references{
Anscombe, F.J. (1948) The transformation of Poisson, binomial and negative-binomial data. \emph{Biometrika} \bold{35(3/4)}, 246--254.

Wang, K., Fingert, J. (2012) Statistical tests for detecting rare variants using variance-stabilizing transformations. Submitted. 
}
\author{
Kai Wang \code{<kai-wang@uiowa.edu>}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{arcsine.test}} and \code{\link{FT.test}}.
}

\examples{
## Example 1 of Li et al. (2010)
G = rbind(c(14, 999), c(3, 1081))
Anscombe.test(G)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
