\name{SNPass.test}
\alias{SNPass.test}
\title{An SNP Association Test using a Proportional Odds Model}
\description{
\code{SNPass.test} performs a test on association between an SNP and one or more phenotypes.
}
\usage{SNPass.test(geno, pheno)}
\arguments{
  \item{geno}{ a length \code{n} vector of SNP genotypes on study subjects. Scoring of the genotypes aa, aA, and AA does not matter as long as the the score for the heterogeneous genotype is between the two homogeneous genotype scores.}
  \item{pheno}{an \code{nxp} matrix of phenotype values. \code{p} could be 1. Phenotype can be dichotomous, continuous, or a mix of both types. }
}

\details{
%% 
The name of this function may change in the future. 
}

\value{
A list with class "\code{htest}" containing the following components:
\item{statistic}{the value of the test statistic.}
\item{p.value}{the p-value for the test computed from a chi-square distribution with 1 df.}
\item{method}{a character string indicting the test performed.}
\item{data.name}{a character string giving the name of the data.}
}


\references{
Wang, K. (2012) Genetic association analysis with multiple phenotypes possibly of different types. Submitted.
}

\author{
Kai Wang \code{<kai-wang@uiowa.edu>}
}

%\note{
%To be added
%}

\examples{
pheno = c(1,1,1,1,1,0,0,0,0,0)
geno = c(0,1,1,2,2,1,1,1,2,2)
SNPass.test(geno, pheno)
}