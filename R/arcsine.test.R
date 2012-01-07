arcsine.test <-
function(G)
# G: (case, ctrl) by (0, 1+2)
{
    n1 = sum(G[1,])
    n2 = sum(G[2,])
    F1 = function(X, n) acos(1-2*X/n)/2    # = asin(2*X/n-1)/2 + pi/2
    S = 4*(F1(G[1,2], n1) - F1(G[2,2], n2))^2 / (1/n1 + 1/n2)
    
	names(S) = "statistic"
    structure(list(statistic = S, p.value = 1-pchisq(S, 1),
                   method = "Comparing Two Proportions Using the Arcsine Transformation", 
                   data.name = deparse(substitute(G))), class = "htest") 
}

