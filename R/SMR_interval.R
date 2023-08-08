#'   Interval Estimates for Summary Data Mendelian Randomization Analysis in the Presence of Winner's Curse
#' 
#'   \code{SMR_interval} calculates conservative box method interval, k-unit support interval, and Wald confidence interval for the causal effect. 
#'   
#' %% ~~ If necessary, more details than the description above ~~ 
#' 
#' @param summary.data a vector (\eqn{\hat b_{gx}}, se(\eqn{\hat b_{gx}}), \eqn{\hat b_{gy}}, se(\eqn{\hat b_{gy}})) of summary data on the exposure \eqn{X} and the outcome \eqn{Y}. Due to winner's curse, the association \eqn{p}-value between the SNP and the exposure is less than \code{sig.level}.
#' @param sig.level the threshold \eqn{p}-value used to select the instrument SNP. The default is \code{5e-8}.
#' @param k the unit used for the k-unit support. The default value is 2.
#' @param alpha \eqn{(1-\alpha)} is the conservative coverage level for the box method interval or the SMR Wald interval. the default value is 0.05 
#' @param method method to construct the interval. It is either "\code{support}", "\code{box}" or "\code{wald}". The default is "\code{box}". 
#' @return The returned value is method-dependent.
#' 
#' For \code{method == "box"}: A list containing the following components:
#' 
#' * an interval estimate.
#' 
#' * type of the interval: completely bounded, exclusive bounded, or bounded.
#' 
#' For \code{method == "support"}: A list containing the following components:
#' 
#' * Estimate The likelihood estimate of \eqn{b}.
#' 
#' * an interval estimate.
#' 
#' For \code{method == "wald"}: an interval estimate.
#' 
#' @author Kai Wang \code{<kai-wang@@uiowa.edu>}
#' @references 
#' Wang, K. (2023) Support interval for two-sample summary data-based mendelia randomization. \emph{Genes}, 14(1):211.
#' 
#' Wang, K. (2023) Interval estimate of causal effect in summary data based Mendelian randomization in the presence of winner’s curse. \emph{Genetic Epidemiology}, 14(1):211.
#' 
#' Zhu, Z. et al. (2016) Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets. \emph{Nature Genetics}, 48(5):481.
#' @examples
#' summary.data = c(0.13707, 0.0235162, -0.0637257, 0.013774)
#' SMR_interval(summary.data)
#' SMR_interval(summary.data, method = "support")
#' SMR_interval(summary.data, method = "wald")
#'
#' @export

SMR_interval = function(summary.data, sig.level = 5e-8, k = 2, alpha = 0.05, method = "box"){
  
  SMR = function(summary.data, alpha = alpha) {         
    x = summary.data[1]/summary.data[2]
    y = summary.data[3]/summary.data[4]

    b = y/x
    V = (1+b^2)/x^2

    c(b-qnorm(1-0.05/2)*sqrt(V), b+qnorm(1-0.05/2)*sqrt(V))*summary.data[4]/summary.data[2]
  }

  box = function(summary.data, alpha = alpha){
    x = summary.data[1]/summary.data[2]
    y = summary.data[3]/summary.data[4]
    tau = qnorm(1-sig.level/2)
    alpha.p = alpha/2   # confidence level for mu_X and mu_Y
  
    box.x = function(x){
      xx = abs(x)

      ff = function(mu.x, sig.level){   
        A = 1-pnorm(tau - mu.x) + pnorm(-tau - mu.x)
        abs((1 - pnorm(xx - mu.x))/A - sig.level)    # Pr(X > x) - sig.level
      }
      
      l = optimize(ff, c(-tau*10, xx - qnorm(1-alpha.p/2)), alpha.p/2)$minimum  
      u = optimize(ff, c(xx-tau*10, xx - qnorm(alpha.p/2)), 1 - alpha.p/2)$minimum   
      
      sort(c(l, u)*sign(x))
    }
  
    box.y = function(y){
      y - c(qnorm(1-alpha.p/2), qnorm(alpha.p/2))
    }
  
    limits.x = box.x(x)
    limits.y = box.y(y)

    if (limits.x[1] <= 0 & limits.x[2] > 0){
      if (limits.y[1] <= 0 & limits.y[2] > 0){
        limits = c(-Inf, Inf)
        class = "unbounded"
      } else if (0 < limits.y[1]) {
        limits = limits.y[1]/limits.x
        class = "exclusive unbounded"
      } else if (limits.y[2] <= 0) {
        limits = limits.y[2]/rev(limits.x)
        class = "exclusive unbounded"
      }
    } else {
      limits = range(outer(limits.y, limits.x, FUN = "/"))
      class = "bounded"
    }

    list(ci = limits*summary.data[4]/summary.data[2], class = class)
  }

  support <- function(summary.data, k=k, sig.level = sig.level){
      lik = function(mu.x, b, x, y, tau) {
        L.x = dnorm(x-mu.x)/(1-pnorm(tau-mu.x) + pnorm(-tau-mu.x))
        L.y = dnorm(y-b*mu.x)
        
        log(L.x*L.y)
      }
      
      l.x = function(mu.x, x, tau) {
        L.x = dnorm(x - mu.x)/(1-pnorm(tau-mu.x) + pnorm(-tau-mu.x))
        log(L.x)
      }
      
      x = summary.data[1]/summary.data[2]
      y = summary.data[3]/summary.data[4]
      s.x = summary.data[2]
      s.y = summary.data[4]
      tau = -qnorm(sig.level/2)
      
      ini.int = c(0, 0)
      if (x > 0) ini.int = c(-1, x)
      if (x < 0) ini.int = c(x, 1)
      mu.hat = optimize(l.x, ini.int, maximum=TRUE, x, tau)$maximum
      b.hat = y/mu.hat
      
      pl = NULL
      if (b.hat >= 0) b.range = seq(y/x-10, b.hat+10, 0.01)
      if (b.hat < 0) b.range = seq(b.hat-10, y/x+10, 0.01)
      for (b0 in b.range){
        ## make sure abs(x-mu.x) < 20 and abs(y-b*mu.x) < 20 to avoid INF in lik function
        if (b0 >= 0) {
          lb = max(x-20, (y-20)/b0)
          ub = min(x+20, (y+20)/b0)
        }
        if (b0 < 0) {
          lb = max(x-20, (y+20)/b0)
          ub = min(x+20, (y-20/b0))
        }
        
        ttt = optimize(lik, c(lb, ub), maximum=TRUE, b0, x, y, tau)
        pl = rbind(pl, c(b0, ttt$maximum, ttt$objective))
      }
      #  plot(pl[,1], pl[,3], xlab = expression(b), ylab = expression(pl(b)))
      
      max.lik = max(pl[,3], na.rm = TRUE)
      indc = which(pl[,3] == max.lik)
      b.hat.2 = pl[indc,1]
      mu.hat.2 = pl[indc, 2]
      support = pl[pl[,3] >= max.lik - k, 1]
      
      list(Estimate = b.hat*s.y/s.x, 
           ci = c(min(support)*s.y/s.x, max(support)*s.y/s.x))  
  }
  
  if (method == "box") return(box(summary.data, alpha = alpha))
  if (method == "support") return(support(summary.data, k=k, sig.level = sig.level))
  if (method == "wald") return(SMR(summary.data, alpha = alpha))
}

