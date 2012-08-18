\name{null.par}
\alias{null.par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the Null Segregation Parameters of the Linear Mixed Model
}
\description{
\code{null.par} estimates the null segregation parameters of the linear mixed model. Its output is used by \code{\link{LMM.test}} to compute the linear mixed model score statistic.
}

\usage{
null.par(K, y, X) 
}

\arguments{
\item{K }{
an \code{nun} BN relatedness matrix computed by function \code{\link{kinf.BN}}. Here \code{n}$ is the number of individuals with non-missing phenotype values.
}
\item{y }{
a numeric vector of length \code{n} containing phenotype values.
}
\item{X }{
an \code{nxp} numeric matrix of covariate values. \code{p} is the number of covariates.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
A profile likelihood is maximized.
}
\value{
A list containing the following components:
\item{y1}{the centered vector of \code{y}.}
\item{Sigma0.1}{the inverse of matrix \code{Sigma0}, the variance matrix of \code{y}.}
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
\code{\link{kinf.BN}} for computing \code{K} and \code{\link{LMM.test}} for conducting the score test based on the linear mixed model.
}

\examples{
G = matrix(sample(c(0,1,2), 200, replace=TRUE), ncol=10)
y = rnorm(10)
X = matrix(rnorm(10), ncol=1)
BN = kinf.BN(G, whole=TRUE)

null.par(BN, y, X)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
