#'   An exact method for SNP-heritability estimation using GWAS summary statistics
#' 
#'   \code{h2_snp} calculates heritability explained by a set of SNPs
#'   
#' %% ~~ If necessary, more details than the description above ~~ 
#' 
#' @param beta a vector of beta coefficients for a set of SNPs. These coefficients are from a GWAS.
#' @param SE a vector of the standard errors of the beta coefficients. 
#' @param N a vector of sample sizes used by the GWAS at these SNPs.
#' @param R LD matrix for these SNPs.
#' @param alpha \eqn{1-\alpha} is the confidence level of the confidence interval.
#' @return A list containing the following components:
#' 
#' * MLE of the heritability.
#' 
#' * umvu (uniformly minimum variance unbiased) estimator of the heritability.
#' 
#' * interval estimate for the heritability.
#' 
#' @author Kai Wang \code{<kai-wang@@uiowa.edu>}
#' @references 
#' Wang, K. (2023) An exact method for SNP-heritability estimation using GWAS summary statistics without heritability modeling. \emph{submitted}
#' @examples
#' beta = c(0.225269, 0.221270, 0.162635, 0.261669, 0.150887, 
#'          0.214515, 0.170296, 0.204454, 0.254811, 0.213803)
#' SE = c(0.033244, 0.032551, 0.032171, 0.031042, 0.032815, 
#'        0.031908, 0.031717, 0.032023, 0.031907, 0.032291)
#' N = 10000
#' R = diag(1, 10)
#' alpha = 0.05
#' h2_snp(beta, SE, N, R, alpha)
#' @export
#' @importFrom MBESS conf.limits.ncf
#' @importFrom MASS ginv
#' @importFrom stats df

h2_snp = function(beta, SE, N, R, alpha){
  f.obj = function(x, TT, df1, df2) df(TT, df1, df2, ncp = x)
  
  mm = nrow(R)  
  nn = round(mean(N),0)
  
  v = beta/sqrt(beta^2 + (N-2)*SE^2)
  ssr.ssy = t(v) %*% ginv(R) %*% v
  TT = ssr.ssy/(1-ssr.ssy)*(nn-1-mm)/mm
  mle = optimize(f.obj, c(0, TT*mm), maximum = TRUE, TT, mm, nn-1-mm)$maximum
  umvu = (nn-3-mm)*ssr.ssy/(1-ssr.ssy) - mm
  ci = conf.limits.ncf(F.value = TT, conf.level = NULL, alpha.lower = alpha/2, alpha.upper = alpha/2, df.1 = mm, df.2 = nn-1-mm)
  
  ttt = c(mle, 
          umvu,
          ci$Lower.Limit, ci$Upper.Limit)
  final = 1/(1+(nn-1)/ttt)
  list(mle = final[1], 
       umvu = final[2], 
       CI = final[3:4])
}



