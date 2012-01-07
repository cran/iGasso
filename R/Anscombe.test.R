Anscombe.test <-
function(G)
# G: (case, ctrl) by (0, 1+2)
{
    n1 = sum(G[1,])
    n2 = sum(G[2,])

    F1 = function(X, n) acos(1-2*X/n)/2    # = asin(2*X/n-1)/2 + pi/2
    F2 = function(X, n)  F1(X+3/8, n+3/4)
    S = 4*(F2(G[1,2], n1) - F2(G[2,2], n2))^2 / (1/(n1+0.5) + 1/(n2+0.5))

	names(S) = "statistic"
    structure(list(statistic = S, p.value = 1-pchisq(S, 1),
                   method = "Comparing Two Proportions Using the Anscombe transformation", 
                   data.name = deparse(substitute(G))), class = "htest") 
}

