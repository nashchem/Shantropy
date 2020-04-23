#' Computes Shannon entropy
#'
#' Computes Shannon entropy for gene expression across cell types or conditions
#' to find cell type or condition specific gene expression patterns
#' genomic regions.
#'
#' Reference papers
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4932457/pdf/FEB4-6-774.pdf
#' https://genomebiology.biomedcentral.com/track/pdf/10.1186/gb-2005-6-4-r33?
#' site=genomebiology.biomedcentral.com
#' https://www.nature.com/nature/journal/v534/n7609/pdf/nature18606.pdf
#'
#' @param data A data frame containing expression data with gene identifier as
#' first column and normalized expression values in RPKM/FPKM, TPM, CPM or any
#' in other columns
#' @param exp_thr A numeric value to filter gene expression
#' @param cell_specificity_thr A numeric value for Shannon entropy threshold to
#' define cellular specificity
#' @return A data frame containing cell type or condition specific expression
#' data in the same format as input data
#' @examples
#' \dontrun{
#' cell_type_data <- compute_ShannonEntropy(data = data, exp_thr = 1,
#' cell_specificity_thr = 1)
#' }
#' @import stats
#' @importFrom dplyr %>%
#' @importFrom assertthat assert_that
#' @export
#'
#'


compute_ShannonEntropy <- function(data = data, exp_thr = 1, cell_specificity_thr = 1){
    # compute total expression
    l <- dim(data)[2]
    exp_rsum <- rowSums(data[, 2:l])

    # computing each gene's individual contribution, Pij
    p_con<- sweep(data[, 2:l], MARGIN = 1, STATS = exp_rsum, FUN = "/")

    # computing shannon's entropy, Hg
    Hg_ind <- (-p_con)*log2(p_con)
    Hg <- rowSums(Hg_ind)

    # Quantify the cellular specificity
    Qg <- Hg - log2(p_con)

    # change colnames
    cell_ids <- colnames(data)[2:l]
    colnames(Qg) <- paste0("Qg_", cell_ids)

    # merge expression and Qg
    zz1 <- cbind(data, Qg)

    # find cell type specific expression data
    # loop over number of cell types
    zz1_cell <- NULL

    for (k in 1:(l-1)){
        zz1_temp <- zz1 %>% dplyr::filter(.[[k+1]] > exp_thr & .[[l+k]] < cell_specificity_thr)
        zz1_temp$cell_id <- cell_ids[k]
        zz1_cell <- rbind(zz1_cell, zz1_temp)
    }

    return(zz1_cell)
}

