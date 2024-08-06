## Statistical test
library(limma)

##################
## Subroutines  ##
##################
limmaTest <- function(data, level, samples, contrasts, comparison, comparisonNames, design, contMatrix) {
    subData = data[, which(colnames(data) %in% samples)]
    nGroups = length(comparison)
    fit = lmFit(data, design) ## Log2-transformation
    fit = contrasts.fit(fit, contMatrix)
    fit = eBayes(fit)
    result = topTable(fit, n = nrow(data), sort = "none", coef = contrasts)

    # Change column names of the result table
    if (nGroups == 2) {
        colnames(result)[which(names(result) == "logFC")] = paste("Log2Fold(", comparisonNames, ")", sep = "")
    } else if (nGroups > 2) {
        ind = grep("AveExpr", colnames(result))
        ind = ind - 1
        for (k in 1:ind) {
            colnames(result)[k] = paste("Log2Fold(", comparisonNames[k], ")", sep = "")
        }
    }
    result$B = NULL
    result$AveExpr = NULL
    result$t = NULL
    result$`F` = NULL
    colnames(result)[which(names(result) == "P.Value")] = "p-value"
    colnames(result)[which(names(result) == "adj.P.Value")] = "FDR"
    colnames(result)[which(names(result) == "adj.P.Val")] = "FDR"
    
    return (list(res = result, data = data))
}

## From isobar package
fitCauchy <- function(x) {
    cauchy.fit <- function(theta, x){
        -sum(dcauchy(x, location = theta[1], scale = theta[2], log = TRUE), na.rm = T)
    }
    good <- !is.na(x) & !is.nan(x)
    x = x[good]
    x = x[x > quantile(x, 0.25) & x < quantile(x, 0.75)]
    theta.start <- c(median(x), IQR(x)/2)
    res <- nlminb(theta.start,cauchy.fit, x = x, lower = c(-10, 1e-20), upper = c(10, 10))
}

cauchyTest <- function(data, level, comparison, comparisonNames) {
    # Assumption: there are only two groups, i.e. two reporters
    subData = data[, which(colnames(data) %in% comparison)]
    log2FC = subData[, 1] - subData[, 2] ## Log2-trasnformed data
    fit = fitCauchy(log2FC)
    pval = sapply(log2FC, function(r) {
        if (is.null(fit) || is.na(log2FC))
            return (NA)
        pcauchy(r, location = fit$par[1], scale = fit$par[2], lower.tail = r < fit$par[1])
    })
    pval = 2 * pval
    fdr = p.adjust(pval, method = "BH")
    result = data.frame(cbind(log2FC, pval, fdr))
    colnames(result) = c(paste0("Log2Fold_", comparisonNames), "p-value", "FDR")

    return (list(res = result, data = subData))
}

statTest = function (data, level, comparison,group,conditions,compareLevel) {
    ## Input arguments
    ##  data: data.frame of id_uni_pep_quan.xlsx or id_uni_prot_quan.xlsx file
    ##  level: analysis level - either "peptide" or "protein"
    ##  comparison: sample names for each group
    ##  group: all sample groups
    ##  condition: conditions that are selected for comparison
    
    ## Retrieve comparison group information
    nGroups = length(comparison)
    
    
    groupinfo <- data.frame(
      row.names = group$Sample,
      condition = factor(conditions)
    )
    
    groups = list()
    nSamples = 0
    samples = NULL
    factors <-levels(groupinfo$condition)
    for (g in 1:nGroups) {
        groups[[g]] = unlist(strsplit(comparison[g], ","))
        nSamples = nSamples + length(groups[[g]])
        samples = c(samples, groups[[g]]) ## Sample labels used in statistical testing
    }
    samples = colnames(data)[colnames(data) %in% samples] ## Sample labels used in statistical testing (in the order of original data)
    
    ## Generate a design matrix (which contains the information of comparison)
    design = model.matrix(~0+groupinfo$condition)
    #design = matrix(0, nrow = nSamples, ncol = nGroups)
    # for (g in 1:nGroups) {
    #     design[which(samples %in% groups[[g]]), g] = 1
    # }
    #colnames(design) = paste("group", seq(1, nGroups), sep = "")
    colnames(design) <- levels(groupinfo$condition)
    
    ## Generate a contrast matrix and new column names for the LIMMA result table
    contVec = NULL
    comparisonNames = NULL
    combMatrix = combn(seq(1, nGroups), 2)
    for (j in 1:ncol(combMatrix)) {
        idx1 = combMatrix[1,j]
        idx2 = combMatrix[2,j]
        contVec = c(contVec, paste(compareLevel[idx1], compareLevel[idx2], sep = "-"))
        comparisonNames = c(comparisonNames, paste(comparison[combMatrix[1, j]], "/", comparison[combMatrix[2, j]], sep = ""))
    }
    contMatrix = makeContrasts(contrasts = contVec, levels = design)
    
    if (nGroups == 2 & max(colSums(design)) == 1) {
        ## Cauchy test
        res = cauchyTest(data, level, comparison, comparisonNames)
    } else if (nGroups > 2 && max(colSums(design)) == 1) {
        stop("For the comparison of multiple groups, replicates are required")
    } else {
        ## LIMMA running
        ## Statistical testing is performed to the "compSamples"
        res = limmaTest(data, level, samples, contVec, comparison, comparisonNames, design, contMatrix)
    }
    return (res)
}
