\name{FT.test}
\alias{FT.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genetic Association Test Based on the Freeman-Tukey Transformation
}
\description{
\code{FT.test} performs a test on association between a variant and case-control status using the Freeman-Tukey transformation.
}
\usage{
FT.test(G)
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
Freeman, M. F., Tukey, J. W. (1950) Transformations related to the angular and the square root. \emph{The Annals of Mathematical Statistics} \bold{21(4)}, 607--611.

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
\code{\link{Anscombe.test}} and \code{\link{arcsine.test}}.
}

\examples{
## Example 1 of Li et al. (2010)
G = rbind(c(14, 999), c(3, 1081))
FT.test(G)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{  }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
