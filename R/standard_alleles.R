#' Determine standard alleles of SNP positions.
#'
#' \code{standard_alleles} determines the reference and variant allele for a SNP position,
#' which are (using total across-sample counts) simply the two most-occurring nucleotides that are also present in the given \code{ref_alleles} column of \code{data_pos}.
#' In case of ties (more than two equally-abundant highest-count alleles or more than one equally-abundant second-highest-count alleles), which should be very rare,
#' the final choice of reference- and/or variant allele is made at random among suitable candidates. 
#'
#' @param data_pos Data frame. Data frame of a SNP position with columns: "chromosome", "position",
#'    "ref_alleles", "dbSNP_ref", "gene", "A", "T", "C", "G", "sample", "sample_nr". At least columns of
#'    allelic counts ("A", "T", "C", "G") and the dbSNP reference alleles ("ref_alleles") should be present.
#'    If no dbSNP reference alleles are available, "A/T/C/G" can be used as reference alleles.
#' @export
#' @return The data as data frame with standard alleles("ref_alleles", "ref" and "var"),
#'     as well as their respective counts ("ref_count" and "var_count").

standard_alleles <- function(data_pos) {
  
  if(is.null(data_pos$ref_alleles[1])){
    data_pos$ref_alleles <- "A/T/C/G"
  }
  
  ref_alleles <- unlist(strsplit(data_pos$ref_alleles[1], "/"))

  allelecount_samples <- lapply(1:nrow(data_pos), function(y) c(data_pos$A[y], data_pos$T[y], data_pos$C[y], data_pos$G[y]))
  allelecount_samples_percent <- lapply(1:nrow(data_pos), function(y) allelecount_samples[[y]] / sum(allelecount_samples[[y]]) * 100)
  allelecount_samples_table <- lapply(1:nrow(data_pos), function(y) setNames(allelecount_samples_percent[[y]], c("A", "T", "C", "G")))

  allelecount_loci <- colSums(matrix(unlist(allelecount_samples), ncol = 4, byrow = T))
  allelecount_loci_percent <- allelecount_loci / sum(allelecount_loci) * 100
  allelecount_loci_table <- setNames(allelecount_loci_percent, c("A", "T", "C", "G"))

  alleles_loci_ordered <- sort(allelecount_loci_table, TRUE)
  if (length(ref_alleles) > 2) {
    standard1 <- NULL
    standard2 <- NULL
    nr <- as.numeric(alleles_loci_ordered)
    char <- names(alleles_loci_ordered)

    standard1 <- which(char %in% ref_alleles)
    standard1 <- char[standard1[which(nr[standard1] == max(nr[standard1]))]]
    if (length(standard1) == 2) {
      data_pos$ref_alleles <- stringr::str_c(standard1, collapse = "/")
    } else if (length(standard1) > 2) {
      data_pos$ref_alleles <- stringr::str_c(sample(standard1, 2), collapse = "/")
    } else {
      nr <- nr[-which(char == standard1)]
      char <- char[-which(char == standard1)]
      standard2 <- which(char %in% ref_alleles)
      standard2 <- sample(char[standard2[which(nr[standard2] == max(nr[standard2]))]], 1)
      ref_alleles <- c(standard1, standard2)
      data_pos$ref_alleles <- stringr::str_c(ref_alleles, collapse = "/")
    }
  }
  st_alleles <- unlist(strsplit(data_pos$ref_alleles[1], "/"))

  data_pos$ref <- names(which.max(allelecount_loci_table[st_alleles]))
  data_pos$var <- st_alleles[-which(st_alleles == data_pos$ref[1])]
  data_pos$ref_count <- data_pos[, which(names(data_pos) == data_pos$ref[1])]
  data_pos$var_count <- data_pos[, which(names(data_pos) == data_pos$var[1])]

  return(data_pos)
}
