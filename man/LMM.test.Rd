\name{LMM.test}
\alias{LMM.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Score Test for the Linear Mixed Model
}
\description{
\code{LMM.test} performs tests on association between an SNP and case-control status. It tests whether the frequencies of an allele are the same between cases and controls. It does not require specification of an inheritance model. 
}

\usage{
LMM.test(g, y1, K, S0.1)
}

\arguments{
\item{g }{
a length \code{n} vector of genotype scores at the SNP being tested. 
}
\item{y1 }{
a numeric vector of length \code{n}. It is the variable \code{y1} returned by function \code{\link{null.par}}. 
}
\item{K }{
an \code{nxn} BN relatedness matrix computed by function \code{\link{kinf.BN}}. Here \code{n} is the number of individuals with non-missing phenotype values.
}
\item{S0.1 }{
an \code{nxn} matrix equal to the inverse of \code{Sigma0}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Matrix multiplication is implemented using existing R functions.
}
\value{
A list with class "\code{test}" containing the following components:
\item{statistic}{the value of the test statistic.}
\item{p.value}{the p-value for the test computed from a 50:50 mixture of 0 and a chi-square distribution with 1 df.}
\item{method}{a character string indicting the test performed.}
\item{data.name}{a character string giving the name of the data.}
}
\references{
Wang K. (2012) A valid and powerful mixed model method for genome-wide association studies. Submitted.
}
\author{
Kai Wang \code{<kai-wang@uiowa.edu>}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{null.par}} for computing \code{y1} and \code{Sigma0.1} and \code{\link{kinf.BN}} for computing \code{K}.
}

\examples{
G = matrix(sample(c(0,1,2), 200, replace=TRUE), ncol=10)
y = rnorm(10)
X = matrix(rnorm(10), ncol=1)
BN = kinf.BN(G, whole=TRUE)
tmp = null.par(BN, y, X)

LMM.test(G[2,], tmp$y1, BN, tmp$Sigma0.1)

apply(G, 1, LMM.test, tmp$y1, BN, tmp$Sigma0.1)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
