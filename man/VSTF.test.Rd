\name{VSTF.test}
\alias{VSTF.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Association Tests for Rare Variants Based on Variance-Stabilizing Transformation
}
\description{
\code{VSTF.test} performs tests on association between a rare variant and case-control status using a variance-stabilizing transformation.
}

\usage{
VSTF.test(G, method = "Anscombe")
}

\arguments{
\item{G }{
a \code{2x2} matrix. The first row is for cases and the second one for controls. In each row, the first element is the number of non-carriers and the second one is the number of carriers with at least 1 copy of the variant. 
}
\item{method }{
a character string indicating which transformation to use. One of \code{"Anscombe"} (default), \code{"arcsine"}, \code{"Freeman-Tukey"}, and \code{"Chanter"}. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Each test is named after the author(s) of the corresponding publication.
}
\value{
A list with class "\code{test}" containing the following components:
\item{statistic}{the value of the test statistic.}
\item{p.value}{the p-value for the test computed from a chi-square distribution with 1 df.}
\item{method}{a character string indicting the test performed.}
\item{data.name}{a character string giving the name of the data.}
}
\references{
Anscombe, F. J. (1948) The transformation of Poisson, binomial and negative-binomial data. \emph{Biometrika} \bold{35(3/4)}, 246--254.

Chanter, D. O. (1975). Modifications of the angular transformation. \emph{Journal of the Royal Statistical Society. Series B (Applied Statistics)} \bold{24(3)}, 354--359.

Freeman, M. F., Tukey, J. W. (1950) Transformations related to the angular and the square root. \emph{The Annals of Mathematical Statistics} \bold{21(4)}, 607--611.

Wang, K., Fingert, J. (2012) Statistical tests for detecting rare variants using variance-stabilizing transformations. \emph{Annals of Human Genetics}. 76(5):402-9.

Zar, J.H. (1999) \emph{Biostatistical Analysis, 4th ed.}, New Jersey:Prentice-Hall, Inc.
}
\author{
Kai Wang \code{<kai-wang@uiowa.edu>}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%%\code{\link{Anscombe.test}} and \code{\link{FT.test}}.
%%}

\examples{
## Example 1 of Li et al. (2010)
G = rbind(c(14, 999), c(3, 1081))
VSTF.test(G)
VSTF.test(G, method = "arcsine")
VSTF.test(G, method = "Freeman-Tukey")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
