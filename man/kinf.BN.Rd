\name{kinf.BN}
\alias{kinf.BN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Balding-Nichols (BN) Relatedness Matrix
}
\description{
\code{kinf.BN} computes a BN relatedness matrix from a matrix of SNP dose scores. 
}
\usage{
kinf.BN(G, whole=FALSE)
}
\arguments{
\item{G }{
a \code{(\# SNP)x(\# individual)} matrix of SNP scores. Each SNP score can be 0, 1, or 2, or imputed dose. 
\code{(\# individual)} includes only individuals with non-missing phenotype values. 
}
\item{whole }{
a logical variable. Are all SNPs included in \code{G}? If not, this function needs to be called more than once. Then the sum of \code{total} is devised by the sum of \code{count}. The default is \code{FALSE}.
}
}
\details{
SNPs with MAF <= 0.01 are excluded. \code{G} can contain missing values. 
}
\value{
If \code{whole = TRUE}, returns the BN relatedness matrix. If \code{whole = FALSE}, returns a list containing the following components:
\item{total}{a matrix equal to the sum of genotype score inner products. Missing genotype scores are excluded.}
\item{count}{a matrix of counts of individuals pairs both genotype scores are non-missing.}
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
\code{\link{null.par}} for computing \code{y1} and \code{Sigma0.1} and \code{\link{LMM.test}} for conducting the score test based on the linear mixed model.
}
\examples{
G1 = matrix(sample(c(0,1,2), 200, replace=TRUE), ncol=10)
G2 = matrix(sample(c(0,1,2), 200, replace=TRUE), ncol=10)

BN = kinf.BN(G1, whole=TRUE)

tmp1 = kinf.BN(G1, whole=FALSE)
tmp2 = kinf.BN(G2, whole=FALSE)
BN = (tmp1$total+tmp2$total)/(tmp1$count+tmp2$count)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
