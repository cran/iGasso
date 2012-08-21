\name{genome.plot}
\alias{genome.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genome-wide Plot of a Variable
}
\description{
\code{genome.plot} plots the value of a variable across the genome.
}

\usage{
genome.plot(mydata, style=1, type="h", sig.line=c(4, -4), 
            sig.color=c("red", "red"), \dots)
}

\arguments{
\item{mydata }{
a data frame containing three variables: \code{y} (numeric, the value of the variable to be plotted), \code{chr} (character, chromosome label), and \code{pos} (numeric, position, for instance, in base pair or centi-Morgan). Examples of \code{y} include -log10 of p-values and test statistic values.
}

\item{style }{
either \code{1} (default) or \code{2}.
}
\item{type }{
a generic graphic parameter. Recommended values are \code{"h"} (default) and \code{"b"}. 
}
\item{sig.line }{
vertical locations of significance lines.
}
\item{sig.color }{
colors of significance lines.
}
\item{\dots }{
other parameters to be passed to function \code{xyplot} in the \code{lattice} package.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function makes use of the function \code{xyplot} from package \code{lattice}. 
}
%%\value{
%%}
%\references{
%}
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
y = rnorm(100)
chr = c(rep(1, 20), rep(3, 20), rep(10, 20), rep(19, 30), rep("X", 10))
pos = c(1:20, 1:20, 1:20, 1:30, 1:10)
mydata = data.frame(y=y, chr=chr, pos=pos)
mydata2 = data.frame(y=y^2, chr=chr, pos=pos)

genome.plot(mydata, sig.line=c(1, -1), ylab="T Statistic")
genome.plot(mydata, sig.line=c(1, -1), ylab="T Statistic", type="b")
genome.plot(mydata2, sig.line=c(2), ylab="y squared")
genome.plot(mydata, style=2, sig.line=c(1, -1), ylab="T Statistic")
genome.plot(mydata, style=2, sig.line=c(1, -1), ylab="T Statistic", type="b")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
