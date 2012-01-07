FT.test <-
function(G)
# G: (case, ctrl) by (0, 1+2)
{
    n1 = sum(G[1,])
    n2 = sum(G[2,])

    F1 = function(X, n) acos(1-2*X/n)/2    # = asin(2*X/n-1)/2 + pi/2
    F3 = function(X, n)  F1(X, n+1) + F1(X+1, n+1)
    S = (F3(G[1,2], n1) - F3(G[2,2], n2))^2 / (1/(n1+0.5) + 1/(n2+0.5))

	names(S) = "statistic"
    structure(list(statistic = S, p.value = 1-pchisq(S, 1),
                   method = "Comparing Two Proportions Using the Freeman-Tukey Transformation", 
                   data.name = deparse(substitute(G))), class = "htest") 
}

