#'   Association Tests for Rare Variants Based on Variance-Stabilizing Transformation
#' 
#'   \code{VSTF.test} performs tests on association between a rare variant and case-control status using a variance-stabilizing transformation.
#'
#' %% ~~ If necessary, more details than the description above ~~ 
#'    Each test is named after the author(s) of the corresponding publication.
#' 
#' @param G a \code{2x2} matrix. The first row is for cases and the second one for controls. In each row, the first element is the number of non-carriers and the second one is the number of carriers with at least 1 copy of the variant. 
#' @param method a character string indicating which transformation to use. One of \code{"Anscombe"} (default), \code{"arcsine"}, \code{"Freeman-Tukey"}, and \code{"Chanter"}.
#' @return A list with class "\code{test}" containing the following components:
#' * statistic the value of the test statistic.
#' * p.value the p-value for the test computed from a chi-square distribution with 1 df.
#' * method a character string indicting the test performed.
#' * data.name a character string giving the name of the data.
#' @author Kai Wang \code{<kai-wang@@uiowa.edu>}
#' @references 
#' Anscombe, F. J. (1948) The transformation of Poisson, binomial and negative-binomial data. \emph{Biometrika} \bold{35(3/4)}, 246--254.
#'
#' Chanter, D. O. (1975). Modifications of the angular transformation. \emph{Journal of the Royal Statistical Society. Series B (Applied Statistics)} \bold{24(3)}, 354--359.
#' 
#' Freeman, M. F., Tukey, J. W. (1950) Transformations related to the angular and the square root. \emph{The Annals of Mathematical Statistics} \bold{21(4)}, 607--611.
#' 
#' Wang, K., Fingert, J. (2012) Statistical tests for detecting rare variants using variance-stabilizing transformations. \emph{Annals of Human Genetics}. 76(5):402-409.
#' 
#' Zar, J.H. (1999) \emph{Biostatistical Analysis, 4th ed.}, New Jersey:Prentice-Hall, Inc.
#' @examples
#' ## Example 1 of Li et al. (2010)
#' G = rbind(c(14, 999), c(3, 1081))
#' VSTF.test(G)
#' VSTF.test(G, method = "arcsine")
#' VSTF.test(G, method = "Freeman-Tukey")
#'
#' @export 

VSTF.test <-
function(G, method="Anscombe")
{
    n1 = sum(G[1,])
    n2 = sum(G[2,])

    F1 = function(X, n) acos(1-2*X/n)/2    # = asin(2*X/n-1)/2 + pi/2
    F2 = function(X, n)  F1(X+3/8, n+3/4)
    F3 = function(X, n)  F1(X, n+1) + F1(X+1, n+1)
    F4 = function(X, n)  { c = 1/(4*n); asin(sqrt(c+(1-2*c)*X/n)) }

    if (method == "arcsine")  S = 4*(F1(G[1,2], n1) - F1(G[2,2], n2))^2/(1/n1 + 1/n2)
    else if (method == "Anscombe")  S = 4*(F2(G[1,2], n1)-F2(G[2,2], n2))^2/(1/(n1+0.5)+1/(n2+0.5))
    else if (method == "Freeman-Tukey")  S = (F3(G[1,2], n1)-F3(G[2,2], n2))^2/(1/(n1+0.5)+1/(n2+0.5))
    else if (method == "Chanter")  S = 4*(F4(G[1,2], n1) - F4(G[2,2], n2))^2/(1/n1 + 1/n2)
    else stop("Valid transformation is arcsine, Anscombe, Freeman-Tukey, or Chanter")

	names(S) = "statistic"
    structure(list(statistic = S, p.value = 1-pchisq(S, 1),
                   method = paste("Comparing Two Proportions Using the", method, "Transformation"), 
                   data.name = deparse(substitute(G))), class = "htest") 
}

