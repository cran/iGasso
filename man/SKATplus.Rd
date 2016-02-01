\name{SKATplus}
\alias{SKATplus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A gene- or pathway-based test of association
}
\description{
\code{SKATlus} provides enhanced power over SKAT by properly estimating the null distribution of SKAT.
}

\usage{
SKATplus(y, G, X=NULL, out_type="D", tau=NULL, permutation=FALSE, B=1000)
}

\arguments{
  \item{y}{a vector of phenotype on \eqn{n} subjects.}
  \item{G}{an \eqn{n\times m} matrix of SNP genotypes of \eqn{n} study subjects at \eqn{m} loci.}
  \item{X}{a matrix of covariates. It has \eqn{n} rows.}
  \item{out_type}{an indicator of the outcome type. \code{"C"} for the continuous outcome and \code{"D"} for the dichotomous outcome.}
  \item{tau}{proportion of selected subjects used for SKATplus.}
  \item{permutation}{an indicator. Use permutation for p-value or not.}
  \item{B}{the number of permutations if \code{permutation=TRUE}}
}


\details{
%%  ~~ If necessary, more details than the description above ~~
This version uses only subjects with lower phenotypic values for estimating the null distribution. That is, the "controls" are those of lower phenotypic values. 
}

\value{
A list with class "\code{htest}" containing the following components:
\item{statistic}{the value of the test statistic, which is the same as SKAT statistic.}
\item{parameter}{sample size and the number of SNPs}
\item{p.value}{the p-value for SKATplus computed using Davies' method.}
\item{method}{a character string indicting the test performed.}
\item{data.name}{a character string giving the name of the data.}
}


\references{
Wang, K. (2016) Boosting the power of the sequence kernel association test (SKAT) almost surely by properly estimating its null distribution. Submitted.
}
\author{
Kai Wang \code{<kai-wang@uiowa.edu>}
}


%\note{
%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%%}

\examples{
n=1000
y = c(rep(1, n/2), rep(0, n/2))
maf = seq(0.05, 0.5, 0.05)
g = NULL
for (j in 1:10){
        geno.freq = c(maf[j]^2, 2*maf[j]*(1-maf[j]), (1-maf[j])^2)
        g = cbind(g, sample(c(0,1,2), n, replace=TRUE, prob=geno.freq))
}
SKATplus(y, g, X=NULL, out_type="D", tau=NULL, permutation=FALSE, B=1000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
