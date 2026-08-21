#' Preprocess Gene Expression Matrix, Filtered low quality genes, and return DESeq normalized matrix
#'
#'
#' @param raw A matrix or data frame containing raw gene expression counts.
#'        Rows represent genes, and columns represent samples. Row names should be 
#'        Ensembl gene IDs.
#' @param validSymbols A data frame containing valid Ensembl gene IDs and HGNC symbols.
#'        Must include the columns `ensembl_gene_id` and `hgnc_symbol`.
#' @param countasZero A numeric threshold for treating counts as zero. Default is 0.
#'        Used for filtering genes based on expression levels.
#'
#' @return A list with the following components:
#'   counts: Matrix of filtered raw counts that passed quality control (QC).
#'   normalized: Matrix of normalized counts after DESeq2 size factor normalization
#'   statistics: statistics for genes and samples

prepMatrix <- function(raw, validSymbols, countasZero = 0, verbose = TRUE){
	stopifnot(!is.null(rownames(raw)))
	stopifnot(is.matrix(raw) | is.data.frame(raw))
	
	filtered <- raw[rowSums(raw > countasZero) > ncol(raw)/3,]
	filtered <- filtered[rownames(filtered) %in% validSymbols$ensembl_gene_id[validSymbols$hgnc_symbol != ""],]
	gsg <- WGCNA::goodSamplesGenes(t(filtered), verbose = verbose)
	filtered <- filtered[gsg$goodGenes,]
	if(!gsg$allOK && verbose){
		if(any(!gsg$goodSamples))
			message("Info: ", paste0(colnames(filtered)[!gsg$goodSamples], collapse = ", "), " did not pass the WGCNA goodSample test. ", 
					"However these samples will not be filtered.")
		if(any(!gsg$goodGenes))
			message("Info: ", sum(!gsg$goodGenes), "Genes did not pass the WGCNA goodSample test after default filtering. ", 
					"However these genes will not be filtered.")
	}
	
	DE <- DESeqDataSetFromMatrix(filtered, colData = data.frame(colnames(filtered)), design = ~1) 
	DE <- estimateSizeFactors(DE)
	normalized <- counts(DE, normalized=TRUE)
    if (verbose)
    	message("#Genes before QC: ",  nrow(raw) , " - #Genes after QC: ", nrow(filtered))
	return(list(counts = filtered, normalized = normalized, 
				statistics = list(genes = data.frame(goodGenes = gsg$goodGenes, 
														   meanCounts = rowMeans(as.matrix(filtered)), 
														   medianCounts = rowMedians(as.matrix(filtered)), 
														   min = apply(filtered, 1,min), max = apply(filtered, 1,max), 
														   symbol = validSymbols$hgnc_symbol[match(rownames(filtered), validSymbols$ensembl_gene_id)])
										,
														   
										samples = data.frame(goodGenes = gsg$goodSamples, 
														   meanCounts = colMeans(as.matrix(filtered)), 
														   medianCounts = colMedians(as.matrix(filtered)), 
														   min = apply(filtered, 2,min), max = apply(filtered, 2,max))
										)
				)
		   )

}

#' Calculate z score of numeric vector
zscore <- function(x)
    return((x - mean(x))/sd(x))

#' Compute mean z-score per sample
#'
#' @param mat A matrix of gene expression values. 
#' @param genes a vector of genes
#'
#' @return vector of mean sample z-scores based on the provided gene module

ZModuleScore <- function(mat, genes){
	inMat_Genes <- genes %in% rownames(mat)
	NA_Genes <- genes[!inMat_Genes]
	if(!is.null(NA_Genes) & length(NA_Genes) > 0) 
		warning("Genes will be excluded because they are not included in matrix: ", 
                paste0(NA_Genes, collapse = ","))
	genes <- genes[inMat_Genes]
	stopifnot(length(genes) > 0)
	mat_subset <- mat[genes,]
	zMean <- t(scale(t(mat_subset))) %>% colMeans
	
	return(zMean)
}


#convert Ensemble <> Symbol
getEns <- function(symb ,conv=symbol_ens){
    inList <- symb %in% symbol_ens$hgnc_symbol
    if (sum(!inList) != 0){
        warning(paste0(symb[!inList], collapse = ","), " not found! Gene will be excluded")
        symb <- symb[inList]
    }
    if (length(symb) == 0)
        return(FALSE)
    
    return(symbol_ens$ensembl_gene_id[match(symb, symbol_ens$hgnc_symbol)])    
}
getSymb <- function(ens ,conv=symbol_ens, ignoreNA = FALSE){
    inList <- ens %in% symbol_ens$ensembl_gene_id
    if (sum(!inList) != 0 & !ignoreNA){
        warning(paste0(ens[!inList], collapse = ","), " not found! Gene will be excluded")
        ens <- ens[inList]
    }
    if (length(ens) == 0 & !ignoreNA)
        return(FALSE)
    
    return(symbol_ens$hgnc_symbol[match(ens, symbol_ens$ensembl_gene_id)])    
}

get_top_TF <- function(data, cond, n = 5){
    if(sum(cond %in% data$condition) == 0)
        stop("provided condition not in data included")
    
    data <- subset(data, condition == cond) %>% arrange(desc(score))
    data <- rbind(head(data), tail(data))
    data <- subset(data, p_value_adj < 0.05)
    return(data)
}