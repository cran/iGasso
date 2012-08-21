\name{MFree.test}
\alias{MFree.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model-free Association Tests
}
\description{
\code{MFree.test} performs tests on association between an SNP and case-control status. It tests whether the frequencies of an allele are the same between cases and controls. It does not require specification of an inheritance model. 
}

\usage{
MFree.test(G, method="score")
}

\arguments{
\item{G }{
a \code{2x3} two-dimensional contingency table in matrix form. The first row is for cases and the second one for controls. In each row, the entries are the number of subjects carrying 0, 1, and 2 copies of the reference allele, respective. 
}
\item{method }{
a character string indicating the test statistic to use. One of \code{"score"} (default), \code{"Wald"}, and \code{"LRT"}. 
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
Wang K. (2012) Statistical tests of genetic association for case-control study designs. \emph{Biostatistics}. Accepted. PMID: 22389176
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
G = rbind(c(161, 474, 489), c(231, 444, 380))
MFree.test(G)
MFree.test(G, method = "Wald")
MFree.test(G, method = "LRT")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
