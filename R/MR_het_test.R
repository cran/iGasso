#'   Test of Heterogeneity in MR using Principal Components
#' 
#'   \code{MR_het_test} performs tests of heterogeneity in MR. 
#'   
#' %% ~~ If necessary, more details than the description above ~~ 
#' 
#' @param x.b a vector of the estimated regression coefficients from the SNP-exposure GWAS 
#' @param y.b a vector of the estimated regression coefficients from the SNP-outcome GWAS
#' @param x.se a vector of SEs for \code{x.b}
#' @param y.se a vector of SEs for \code{y.b}
#' @param b0 a value used for the common effect size. It is used for the weighting matrix
#' @param k the number of principal components used. It is used by the \eqn{\tilde Q(b_0)} statistic. The default is NULL
#' @param cum.prop threshold for selecting \code{k}. It is void if \code{k} is specified. The default is 0.8, i.e., the proportion of variance explained by the top \code{k} principal components is 0.8
#' @return A list containing the following components:
#' 
#' * \eqn{P_\text{min}(b_0)} statistic and its \eqn{P}-value.
#' 
#' * \eqn{\tilde Q_\text{min}(b_0)} statistic, its degrees of freedom, and its \eqn{P}-value.
#' 
#' @author Kai Wang \code{<kai-wang@@uiowa.edu>}
#' @references 
#' Wang, K, Alberding, Steven Y. (2024) Powerful test of heterogeneity in two-sample summary-data Mendelian randomization. Submitted.
#'
#' @examples
#' p = 10
#' b = 0.5
#' gamma.true = runif(p, 0.34, 1.1)
#' x.se = runif(p, 0.06, 0.1)
#' y.se = runif(p, 0.015, 0.11)
#' x.b = rnorm(p, gamma.true, x.se)
#' y.b = rnorm(p, b*gamma.true, y.se)
#' b0 = 0.4
#' 
#' MR_het_test(x.b, y.b, x.se, y.se, b0)
#'
#' @export

MR_het_test = function(x.b, y.b, x.se, y.se, b0, k = NULL, cum.prop = 0.8){
  p = length(x.b)
  Sigma = diag(y.se^2 + b0^2*x.se^2)  
  Sigma.1 = solve(Sigma)
  x.b.Sigma.1 = as.vector(crossprod(x.b, Sigma.1))
  dot.Sigma.inv = Sigma.1 - outer(x.b.Sigma.1, x.b.Sigma.1) / sum(x.b * x.b.Sigma.1)
  pca = eigen(dot.Sigma.inv)
  A = pca$vectors[,-p]

  z.stat.2 = crossprod(A, y.b - b0*x.b)^2*pca$values[-p]
  max.z.stat.2 = max(z.stat.2)
  min.p = 1 - pchisq(max.z.stat.2, 1)
  P.min.stat = data.frame(P.min = min.p, p.value = 1 - (1 - min.p)^(p-1))

  if (is.null(k)) {
    k = min((1:ncol(Sigma))[cumsum(pca$values)/sum(pca$values) > cum.prop])
  }
  Q.tilde = sum(((y.b - b0*x.b) %*% pca$vectors[,1:k])^2 * pca$values[1:k])
  Q.tilde.stat = data.frame(Q.tilde.min = Q.tilde, k = k, p.value = 1-pchisq(Q.tilde, k)) 

  list(P.min.test = P.min.stat, Q.tilde.min.test = Q.tilde.stat)  
}

