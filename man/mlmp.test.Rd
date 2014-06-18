\name{mlmp.test}
\alias{mlmp.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A suite of multilocus multiple phenotype tests for association
}
\description{
\code{mlmp.test} performs tests of association between multilocus SNPs and one or more traits possibly of different types.
}

\usage{
mlmp.test(g, y, weights = NULL, stat = "score")
}

\arguments{
  \item{g}{an \eqn{n\times p} matrix of SNP genotypes of \eqn{n} study subjects at \eqn{p} loci.}
  \item{y}{an \eqn{n\times q} matrix of phenotype values of \eqn{n} study subjects on \eqn{q} traits. \eqn{q} could be 1. Traits can be dichotomous, continuous, selected data, etc. or a mix of these. }
  \item{weights}{a vector of length \eqn{p}. Each element is the weight used for the corresponding variant. The default weight for each variant is the inverse of the sample variance of the genotype score at this variant. This option differs with the SKAT option by the same name in that its elements are \eqn{w}, not the square root of \eqn{w}. }
  \item{stat}{Test statistic to be used. One of \code{"F"}, \code{"Wald"}, or \code{"score"}(default).}
}


\details{
%%  ~~ If necessary, more details than the description above ~~
This method regresses multilocus genotype over multiple phenotypes and test a quadratic null hypothesis. The \eqn{p}-value of a test statistic is determined through a linear combination of independent chi-square distributions and is evaluated via Davies' method implemented in package \code{CompQuadForm}.
}

\value{
A list with class "\code{htest}" containing the following components:
\item{statistic}{the value of the test statistic.}
\item{parameter}{the number of SNPs and the number of traits.}
\item{p.value}{the p-value for the test computed using Davies' method.}
\item{method}{a character string indicting the test performed.}
\item{data.name}{a character string giving the name of the data.}
}


\references{
Wang, K. (2014) Multilocus genetic association analysis with multiple phenotypes. Submitted.
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
n=400
y = c(rep(1, n/2), rep(0, n/2))
y = cbind(y, rnorm(n))
y = cbind(y, rnorm(n))

maf = seq(0.05, 0.5, 0.05)
g = NULL
for (j in 1:10){
        geno.freq = c(maf[j]^2, 2*maf[j]*(1-maf[j]), (1-maf[j])^2)
        g = cbind(g, sample(c(0,1,2), n, replace=TRUE, prob=geno.freq))
}
mlmp.test(g, y, weights=rep(1,10))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
