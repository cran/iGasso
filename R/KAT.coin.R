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
