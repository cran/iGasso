#' Conditional Inference for the Kernel Association Test (KAT)
#' 
#' Computes the asymptotic and the approximate conditional p-values for the kernel association test
#' 
#' %% ~~ If necessary, more details than the description above ~~ 
#'     The asymptotic conditional null distribution is obtained using results in Strasser and Weber (1999). 
#'     The p-value based on this distribution is computed using Davies' method.
#' 
#' @param y a vector of phenotype on \eqn{n} subjects.
#' @param G an \eqn{n \times m} matrix of SNP genotypes of \eqn{n} study subjects at \eqn{m} loci.
#' @param X a matrix of covariates. It has \eqn{n} rows.
#' @param out_type an indicator of the outcome type. \code{"C"} for the continuous outcome and \code{"D"} for the dichotomous outcome.
#' @param distribution a character, the conditional null distribution of the test statistic can be approximated by its asymptotic distribution ("asymptotic", default) or via Monte Carlo resampling ("approximate"), as in package \code{coin}.
#' @param B the number of permutations if \code{distribution = "approximate"}.
#' @return A list with class "\code{htest}" containing the following components:
#' 
#' * statistic the value of the kernel association test statistic.
#' 
#' * parameter sample size and the number of SNPs.
#' 
#' * p.value the p-value based on the asymptotic or the approximate conditional null distribution.
#' 
#' * method a character string indicting the test performed.
#' 
#' * data.name a character string giving the name of the data.
#' @author Kai Wang \code{<kai-wang@@uiowa.edu>}
#' @references 
#' Strasser, H. and Weber, C. (1999) On the asymptotic theory of permutation statistics. \emph{Mathematical Methods of Statistics}. 8(2):220-250.
#' 
#' Wang, K. (2017) Conditional Inference for the Kernel Association Test. Bioinformatics 33 (23), 3733-3739.
#' @keywords conditional inference
#' @examples
#' n=1000
#' y = c(rep(1, n/2), rep(0, n/2))
#' maf = seq(0.05, 0.5, 0.05)
#' g = NULL
#' for (j in 1:10){
#'    geno.freq = c(maf[j]^2, 2*maf[j]*(1-maf[j]), (1-maf[j])^2)
#'    g = cbind(g, sample(c(0,1,2), n, replace=TRUE, prob=geno.freq))
#'    }
#' KAT.coin(y, g, X=NULL, out_type="D", B=1000)
#'
#' @export 
#' @importFrom CompQuadForm davies

KAT.coin = function (y, G, X = NULL, out_type = "D", distribution = "asymptotic", B = 1000) 
{
    DNAME = paste(deparse(substitute(y)), "and", deparse(substitute(G)), "and", deparse(substitute(X)))
    if (deparse(substitute(X))!="NULL") DNAME = paste(DNAME, "and", deparse(substitute(X))) 
    n = length(y)
    G = as.matrix(G, nrow = n)
    if (is.null(X)) {
    	    X = matrix(1, ncol = 1, nrow = n)
    	    m = 0
    	}
    else {
    	    X = cbind(1, X)
    	    m = ncol(X)
    	}
    k = ncol(X)
    if (out_type == "D") fit0 = glm(y ~ X - 1, family = "binomial")
    if (out_type == "C") fit0 = lm(y ~ X - 1)
    res.y = y-fitted(fit0)

    stat= sum(crossprod(G, res.y)^2)
    if (distribution == "asymptotic"){
#    	    pca0 = (n-1)*drop(var(res.y))*svd(cov(G), nu=0, nv=0)$d
    	    pca0 = (n-m-1)*drop(var(res.y))*svd(cov(G), nu=0, nv=0)$d
        p.value = davies(1, lambda = pca0/stat)$Qq
        method = "SKAT test of association with conditional asymptotic p-value"
    }
    if (distribution == "approximate"){
	    stat.p = 0
        for (i in 1:B) stat.p = stat.p + (sum(crossprod(G, sample(res.y))^2) >=stat)
        p.value = stat.p/B
        method = "SKAT test of association with conditional approximate p-value"
    }

    PAR = c(n, ncol(G))
    names(PAR) = c("#subjects", "#SNPs")
    names(stat) = "SKAT"
    structure(list(statistic = stat, p.value = p.value, parameter = PAR, method = method, 
        data.name = DNAME), class = "htest")
}	
