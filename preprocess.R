preprocess = function(df, level, metric, pct) {
    # Input arguments
    # df = dataframe from either "id_uni_prot_quan.xlxs" or "id_uni_pep_quan.xlsx"
    # level = either "peptide" or "protein"
    # metric = CV or MAD (for the selection of highly variable features)
    # pct = percentage of the above metric (e.g., 10 = top 10% of features according to the above metric)
    
    if (level == "peptide") {
        entry = df[, 1]
    } else if (level == "protein") {
        entry = row.names(df)
    }

    # colInd = grep("sig", colnames(df))
    colInd = c(3:ncol(df))
    exprs = log(df[, colInd], 2)
    if (!is.null(metric) & !is.null(pct)) {
        # For exploratory analysis
        # CV or MAD calculation is based on log2-transformed intensities, but output format is raw-intensity scale
        cv = apply(exprs, 1, sd, na.rm = T) / rowMeans(exprs, na.rm = T)
        mad = apply(abs(exprs - apply(exprs, 1, median, na.rm = T)), 1, median, na.rm = T)
        threshold = as.numeric(pct)/ 100 ## Threshold percentage
        rowInd = NULL
        if (as.numeric(metric) == 1) {
            rowInd = cv > quantile(cv, prob = 1 - threshold, na.rm = T)
        } else if (as.numeric(metric) == 2) {
            rowInd = mad > quantile(mad, prob = 1 - threshold, na.rm = T)
        }
    } else {
        # For DE analysis
        rowInd = c(1: nrow(exprs))
    }
    
    # Organize a dataset for subsequent analyses
    # Column 1: entry (either protein accession or peptide sequence)
    # Column 2~ : log2-transformed intensity values for reporter ions
    entry = entry[rowInd]
    res = log(df[rowInd, colInd], 2)
    rownames(res) = entry
    
    # Organize a dataset for downloading
    # colInd = max(grep("sig", colnames(df)))
    # df = df[, 1:colInd] ## Remove statistical analysis results from previous jump -q
    # df = df[rowInd, ]
    return (list(rawData = df, data = res, level = level))
}