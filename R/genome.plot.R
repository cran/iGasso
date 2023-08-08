#'   Genome-wide Plot of a Variable
#' 
#'   \code{genome.plot} plots the value of a variable across the genome.
#'
#' %% ~~ If necessary, more details than the description above ~~ 
#'        This function makes use of the function \code{xyplot} from package \code{lattice}. 
#' 
#' @param mydata a data frame containing three variables: \code{y} (numeric, the value of the variable to be plotted), \code{chr} (character, chromosome label), and \code{pos} (numeric, position, for instance, in base pair or centi-Morgan). Examples of \code{y} include -log10 of p-values and test statistic values.
#' @param style either \code{1} (default) or \code{2}.
#' @param type a generic graphic parameter. Recommended values are \code{"h"} (default) and \code{"b"}. 
#' @param sig.line vertical locations of significance lines.
#' @param sig.color colors of significance lines.
#' @param \dots other parameters to be passed to function \code{xyplot} in the \code{lattice} package.
#' @author Kai Wang \code{<kai-wang@@uiowa.edu>}
#' @examples
#'   y = rnorm(100)
#'   chr = c(rep(1, 20), rep(3, 20), rep(10, 20), rep(19, 30), rep("X", 10))
#'   pos = c(1:20, 1:20, 1:20, 1:30, 1:10)
#'   mydata = data.frame(y=y, chr=chr, pos=pos)
#'   mydata2 = data.frame(y=y^2, chr=chr, pos=pos)
#'   
#'   genome.plot(mydata, sig.line=c(1, -1), ylab="T Statistic")
#'   genome.plot(mydata, sig.line=c(1, -1), ylab="T Statistic", type="b")
#'   genome.plot(mydata2, sig.line=c(2), ylab="y squared")
#'   genome.plot(mydata, style=2, sig.line=c(1, -1), ylab="T Statistic")
#'   genome.plot(mydata, style=2, sig.line=c(1, -1), ylab="T Statistic", type="b")
#'
#' @export 
#' @importFrom lattice xyplot

genome.plot = function(mydata, style=1, type="h", sig.line=c(4, -4), sig.color=c("red", "red"), ...)
{
	mydata = na.omit(mydata)
    t.chr = as.character(mydata$chr)
    t.chr = ifelse(t.chr %in% c("1","2","3","4","5","6","7","8","9"), paste("0", t.chr, sep=""), t.chr)
    t.chr = factor(t.chr)
    chr.name = levels(t.chr)
    chr.name = ifelse(substr(chr.name, 1,1) == "0", substring(chr.name, 2), chr.name)
    levels(t.chr) = chr.name
    n.chr = length(chr.name)
    mins = as.vector(tapply(mydata$pos, t.chr, FUN="min"))
    maxs = as.vector(tapply(mydata$pos, t.chr, FUN="max"))

    if (style == 1){
    	xyplot(y ~ pos|t.chr, data=mydata, xlab = "Chromosome", type=type, ..., 
               layout=c(n.chr,1), scales=list(x=list(relation="free", draw=FALSE)), 
               par.settings = list(layout.widths=list(panel=maxs-mins), axis.line=list(lwd=0.1),                
                                   strip.border=list(lwd=0.1)),
               strip = function(..., bg, par.strip.text) 
                               strip.default(..., bg="pink",  par.strip.text=list(cex=0.75)),
               abline=list(h=sig.line, col=sig.color))
    }
    else if(style == 2){
    	xyplot(y ~ pos|t.chr, data=mydata, xlab = "Chromosome", type=type, ..., 
               layout=c(n.chr,1), strip=FALSE,
               scales=list(x=list(relation="free", tck=c(0,0),
                           at=as.vector((maxs+mins)/2, mode="list"), 
                           labels=as.vector(chr.name, mode="list"))), 
               par.settings = list(layout.widths=list(panel=maxs-mins), axis.line=list(lwd=0.1)), 
               abline=list(h=sig.line, col=sig.color))
    }
    else stop("The value of shape should be either 1 or 2")
}



