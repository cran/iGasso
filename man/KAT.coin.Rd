\name{KAT.coin}
\alias{KAT.coin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conditional Inference for the Kernel Association Test (KAT)
}
\description{
\code{KAT.coin} computes the asymptotic and the approximate conditional p-values for the kernel association test
}

\usage{
KAT.coin(y, G, X = NULL, out_type = "D", distribution = "asymptotic", B = 1000) 
}

\arguments{
  \item{y}{a vector of phenotype on \eqn{n} subjects.}
  \item{G}{an \eqn{n\times m} matrix of SNP genotypes of \eqn{n} study subjects at \eqn{m} loci.}
  \item{X}{a matrix of covariates. It has \eqn{n} rows.}
  \item{out_type}{an indicator of the outcome type. \code{"C"} for the continuous outcome and \code{"D"} for the dichotomous outcome.}
  \item{distribution}{a character, the conditional null distribution of the test statistic can be approximated by its asymptotic distribution ("asymptotic", default) or via Monte Carlo resampling ("approximate"), as in package \code{coin}.}
  \item{B}{the number of permutations if \code{distribution = "approximate"}}
}


\details{
%%  ~~ If necessary, more details than the description above ~~
The asymptotic conditional null distribution is obtained using results in Strasser and Weber (1999). The p-value based on this distribution is computed using Davies' method.
}

\value{
A list with class "\code{htest}" containing the following components:
\item{statistic}{the value of the kernel association test statistic.}
\item{parameter}{sample size and the number of SNPs}
\item{p.value}{the p-value based on the asymptotic or the approximate conditional null distribution.}
\item{method}{a character string indicting the test performed.}
\item{data.name}{a character string giving the name of the data.}
}


\references{
Strasser, H. and Weber, C. (1999) On the asymptotic theory of permutation statistics. \emph{Mathematical Methods of Statistics}. 8(2):220-250.

Wang, K. (2016) Conditional Inference for the Kernel Association Test. Submitted.
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
KAT.coin(y, g, X=NULL, out_type="D", B=1000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
