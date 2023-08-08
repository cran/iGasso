#'   A Gene- or Pathway-Based Test of Association
#'   
#'   \code{SKATlus} provides enhanced power over SKAT by properly estimating the null distribution of SKAT.
#'
#' %% ~~ If necessary, more details than the description above ~~ 
#'    This version uses only subjects with lower phenotypic values for estimating the null distribution. That is, the "controls" are those of lower phenotypic values. When "controls" are of higher phenotypic values, change the sign of the phenotypic values in order to use this function.
#' 
#' @param y a vector of phenotype on \eqn{n} subjects.
#' @param G an \eqn{n \times m} matrix of SNP genotypes of \eqn{n} study subjects at \eqn{m} loci.
#' @param X a matrix of covariates. It has \eqn{n} rows.
#' @param out_type an indicator of the outcome type. \code{"C"} for the continuous outcome and \code{"D"} for the dichotomous outcome.
#' @param tau proportion of selected subjects used for SKATplus.
#' @param permutation an indicator. Use permutation for p-value or not.
#' @param B the number of permutations if \code{permutation=TRUE}
#' 
#' @return   A list with class "\code{htest}" containing the following components:
#' * statistic the value of the test statistic, which is the same as SKAT statistic.
#' * parameter sample size and the number of SNPs
#' * p.value the p-value for SKATplus computed using Davies' method.
#' * method a character string indicting the test performed.
#' * data.name a character string giving the name of the data.
#' @author Kai Wang \code{<kai-wang@@uiowa.edu>}
#' @references 
#'     Wang, K. (2016) Boosting the power of the sequence kernel association test (SKAT) almost surely by properly estimating its null distribution. \emph{A J Hum Genet}. 99 (1), 104-114.
#'
#' @examples
#'     n=1000
#'     y = c(rep(1, n/2), rep(0, n/2))
#'     maf = seq(0.05, 0.5, 0.05)
#'     g = NULL
#'     for (j in 1:10){
#'        geno.freq = c(maf[j]^2, 2*maf[j]*(1-maf[j]), (1-maf[j])^2)
#'        g = cbind(g, sample(c(0,1,2), n, replace=TRUE, prob=geno.freq))
#'        }
#'        SKATplus(y, g, X=NULL, out_type="D", tau=NULL, permutation=FALSE, B=1000)
#'
#' @export 

SKATplus = function(y, G, X=NULL, out_type="D", tau=NULL, permutation=FALSE, B=1000){
	DNAME = paste(deparse(substitute(G)), "and", deparse(substitute(y)), "and", deparse(substitute(X)))
    n = length(y)
    G = as.matrix(G, nrow = n)      # genotype matrix
	if(is.null(X)) X = matrix(1, ncol=1, nrow=n) else X = cbind(1, X)
    k = ncol(X)                               # k=p+1

    if (out_type == "D"){
     	fit0 = glm(y ~ X, family="binomial")
    	    hpi = fitted.values(fit0)
    	    res.y = y - hpi
    	    V  = hpi*(1-hpi)
    }
    if (out_type == "C"){
       	fit0 = lm(y ~ X)
	    res.y = residuals(fit0)
	    	V = rep(summary(fit0)$sigma^2, n)
    }
    stat = sum(crossprod(G, res.y)^2)

    if (is.null(tau)) {
    	    selected = (res.y <= 0) 
    	    tau = mean(selected)
    	    }
    	else selected = (rank(res.y) <= tau*n)
    G0 = G[selected,]  
    n0 = sum(selected)
    if (permutation){
        stat.p = rep(1,B)
        for (i in 1:B) {
          	res.y0 = sample(res.y, n0)
          	stat.p[i] = sum(crossprod(G0, res.y0 - mean(res.y0))^2)*(n-k)/(n0-k)
        }
        pvalue = mean(stat.p>=stat)
    }
    else{
        
        PPP = function(V, X){
     	    XV = crossprod(X, diag(V))
    	    diag(V) - crossprod(XV, solve(XV %*% X, XV))
        }      

     	P0 = PPP(V[selected], X[selected,])
    	    GPG0 = crossprod(G0, P0 %*% G0)
        pca0 = eigen(GPG0, symmetric=TRUE)
        pvalue = davies(stat, lambda=pca0$values*(n-k)/(n0-k))$Qq
    }

    PAR = c(ncol(G), n, round(tau*100,0))
    names(PAR) = c("#SNPs", "#subjects", "% used for the null")
    names(stat) = "SKAT"
    structure(list(statistic = stat, p.value = pvalue, parameter = PAR, 
        method = "SKAT+: gene- or pathway-based test of association based on a proper null", data.name = DNAME), 
        class = "htest")
}
