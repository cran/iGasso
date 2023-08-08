#' An Accurate Normalization Method for High Throughput Sequencing Data
#' 
#' Estimates scaling factors using the trimmed average of ratios of quantiles (TARQ) method
#' 
#' %% ~~ If necessary, more details than the description above ~~ 
#' Estimation of scaling factors for NGS read counts data is challenging. TARQ provides a quantile-based method for estimating 
#' scaling factors. It starts by ordering the raw counts 
#' sample by sample and constructs a reference sample from these ordered counts. To compute the scaling factor for a sample, ratios
#' of its quantiles to those of the reference sample are formed. Zero ratios are removed. Then extreme ratios (too large or too 
#' small) are trimmed before taking average over the remaining ratios. 
#' 
#' @param X a matrix of raw counts. Rows are for taxa (genes, transcripts) and columns for samples
#' @param tau a numerical value in (0, 0.5). The upper \eqn{\tau/2 \times 100}\% and the lower \eqn{\tau/2 \times 100}\% of the 
#' ratios of quantiles are trimmed 
#' @return a vector of scaling factors. Normalized counts can be obtained by \code{sweep(X, 2, scale.factors, FUN="/")}
#' @author Kai Wang \code{<kai-wang@@uiowa.edu>}
#' @references Wang, K. (2018) An Accurate Normalization Method for Next-Generation Sequencing Data.
#' Submitted.
#' @keywords Next-Generation Sequencing RNA-Seq microbiome normalization
#' @examples
#' 
#' #data(throat.otu.tab)
#' #data(throat.meta)
#' #otu.tab = t(throat.otu.tab)
#' #tarq(otu.tab, 0.3)
#' 
#' ##### Use TARQ with DESeq2 
#' #dds <- DESeqDataSetFromMatrix(countData = otu.tab,
#'  #                             colData = throat.meta,
#'  #                             design= ~ SmokingStatus)
#' #sizeFactors(dds) <- tarq(otu.tab, 0.3)
#' #dds <- DESeq(dds)                 
#' #results(dds)
#' #
#' ###### Use TARQ with edgeR
#' #cs <- colSums(otu.tab)
#' #scale.factors <- tarq(otu.tab, 0.3)
#' #tmp <- scale.factors/cs
#' #norm.factors <- tmp/exp(mean(log(tmp)))
#' #dgList <- DGEList(counts = otu.tab, genes=rownames(otu.tab), norm.factors = norm.factors)
#' #designMat <- model.matrix(~ throat.meta$SmokingStatus)
#' #dgList <- estimateGLMCommonDisp(dgList, design=designMat)
#' #fit <- glmFit(dgList, designMat)
#' #glmLRT(fit, coef=2)
#'
#' @export tarq
#' @importFrom stats quantile

tarq = function(X, tau=0.3){
    X = apply(X, 2, FUN="sort")
#    ref = rowSums(X)
    ref = exp(rowMeans(log(X), na.rm=TRUE))
	  alpha = function(y) { 
    		idx = (1:length(y))[y > 0]               # keep non-zero counts
    		r = y[idx]/ref[idx]                          # ratio of quantiles
		    r = r[r >= quantile(r, tau/2) & r <= quantile(r, 1-tau/2)]
        return(exp(mean(log(r)))) 
	  }
    
    scale.factor =  apply(X, 2, FUN=alpha)
    scale.factor
}
