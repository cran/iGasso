#' Model-free Association Tests
#' 
#' \code{MFree.test} performs tests on association between an SNP and case-control status. It tests whether the frequencies of an allele are the same between cases and controls. It does not require specification of an inheritance model. 
#'
#' %% ~~ If necessary, more details than the description above ~~ 
#'    Each test is named after the author(s) of the corresponding publication.
#' 
#' @param G a \code{2x3} two-dimensional contingency table in matrix form. The first row is for cases and the second one for controls. In each row, the entries are the number of subjects carrying 0, 1, and 2 copies of the reference allele, respective. 
#' @param method a character string indicating the test statistic to use. One of \code{"score"} (default), \code{"Wald"}, and \code{"LRT"}. 
#' @return   A list with class "\code{test}" containing the following components:
#' * statistic the value of the test statistic.
#' * p.value the p-value for the test computed from a chi-square distribution with 1 df.
#' * method a character string indicting the test performed.
#' * data.name a character string giving the name of the data.
#' @author Kai Wang \code{<kai-wang@@uiowa.edu>}
#' @references 
#' Wang K. (2012) Statistical tests of genetic association for case-control study designs. \emph{Biostatistics}. 13(4):724-33. PMID: 22389176
#' @examples
#' G = rbind(c(161, 474, 489), c(231, 444, 380))
#' MFree.test(G)
#' MFree.test(G, method = "Wald")
#' MFree.test(G, method = "LRT")
#'
#' @export 

MFree.test <-
function(G, method="score")
{
    n1 = G[1,]
    n2 = G[2,]


	l1 = function(n1, n2)
	sum(n1[n1>0]*log(n1[n1>0]/sum(n1))) + sum(n2[n2>0]*log(n2[n2>0]/sum(n2)))

    l0 = function(p, n1, n2)
    {
	    q = 1-p
	    gamma1 = F.est(p, n1)
	    gamma2 = F.est(p, n2)
	    P1 = c(q^2 + gamma1, 2*p*q - 2*gamma1, p^2 + gamma1)
	    P2 = c(q^2 + gamma2, 2*p*q - 2*gamma2, p^2 + gamma2)

    	sum(n1[n1>0]*log(P1[n1>0])) + sum(n2[n2>0]*log(P2[n2>0]))
    }

    LR.test = function(n1, n2)
    {
	    ll0 = optimize(l0, c(0, 1), n1=n1, n2=n2,tol = 1e-04, maximum=TRUE)
	    2*l1(n1, n2)-2*ll0$objective
    }

    pf = function(n)
    {
	    p = n[3]/sum(n) + n[2]/sum(n)/2
        gamma = p*(1-p)-n[2]/sum(n)/2	   # gamma = F*p*(1-p)
        c(p, gamma)
    }

    F.est = function(p, n)
    # maximum likelihood estimation of gamma := F*p*(1-p) for given value of freq p
    {
	    q = 1-p
        shift = 0.0000000001

    	lb = max(c(-p^2, -q^2, p*q-0.5))+shift
	    ub = min(c(1-p^2, 1-q^2, p*q))-shift

    	if (n[1]>0 & n[2]>0 & n[3]>0)
	    {
		    a = sum(n)
		    b = (n[1]+n[2])*p^2+(n[3]+n[2])*q^2-(n[1]+n[3])*p*q
		    cc = (n[2]*p*q - n[1]*p^2 - n[3]*q^2)*p*q
		    temp = (-b+sqrt(b^2-4*a*cc))/2/a   # desired root of the first-order equation
	    }
	    else if (n[1]<1 & n[2]>0 & n[3]>0) temp = n[3]*p/sum(n) - p^2
	    else if (n[1]>0 & n[2]>0 & n[3]<1) temp = n[1]*q/sum(n) - q^2
	    else if (n[1]>0 & n[2]<1 & n[3]>0) temp = ub
	    else if (n[1]>0 & n[2]<1 & n[3]<1) temp = ub
	    else if (n[1]<1 & n[2]>0 & n[3]<1) temp = lb
	    else if (n[1]<1 & n[2]<1 & n[3]>0) temp = ub

        temp = max(temp, lb)
	    temp = min(temp, ub)
       	est = c(lb, ub, temp)
        obj = 0
        if (n[1]>0) obj = obj +n[1]*log(q^2+est)
        if (n[2]>0) obj = obj +n[2]*log(2*p*q-2*est)
        if (n[3]>0) obj = obj +n[3]*log(p^2+est)

    	est[order(obj)][length(est)]   
    }

    part = function(n, p)
    {
	    gamma = F.est(p, n)
	    q = 1-p
	    P.tilde = c(q^2+gamma, 2*p*q-2*gamma, p^2+gamma)
	    P.hat = n/sum(n)
    		
	    sum((P.hat/P.tilde - 1)*n)
    }

    score.test = function(n1, n2)
    {
	    p = pf(n1+n2)[1]

        part(n1, p)+part(n2, p)
    }

    Wald.test = function(n1, n2)
    {
    	p1 = pf(n1)[1]
	    gamma1 = pf(n1)[2]            # gamma1 = F1*p1*(1-p1)
		p2 = pf(n2)[1]
		gamma2 = pf(n2)[2]            # gamma1 = F2*p2*(1-p2)
    	p = pf(n1+n2)[1]
    	gamma = pf(n1+n2)[2]          # gamma1 = F*p*(1-p)
    	SigmaW = (p1*(1-p1)+gamma1)/sum(n1) + (p2*(1-p2)+gamma2)/sum(n2)
    	
    	2*(p1-p2)^2/SigmaW
    }


    if (method == "score")  S = score.test(n1, n2)
    else if (method == "Wald")  S = Wald.test(n1, n2)
    else if (method == "LRT")  S = LR.test(n1, n2)
    else stop("Valid test is score, Wald, or LRT")

	names(S) = "statistic"
    structure(list(statistic = S, p.value = 1-pchisq(S, 1),
                   method = paste("Model-Free test of Association:", method, "Statistic"), 
                   data.name = deparse(substitute(G))), class = "htest") 
}



