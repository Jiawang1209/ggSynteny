#' Convert GFF3 annotation to MCScanX gff format
#'
#' @param gffpath Character
#' Path to the input GFF3 annotation file.
#' @param output Character
#' Path to the output file in MCScanX gff format.
#'
#' @returns None
#' @export
#'
#' @examples NULL
parse_gff3_to_mcscanx <- function(gffpath, output){
  
  gff <- readr::read_delim(file = gffpath, delim = "\t", comment = "#", col_names = F) %>%
    dplyr::filter(X3 == "gene") %>%
    dplyr::select(X1, X4, X5, X9) %>%
    dplyr::mutate(X9 = str_remove(X9, pattern = ";.*")) %>%
    dplyr::mutate(X9 = str_remove(X9, pattern = "ID=")) %>%
    dplyr::select(1,4,2,3)
  
  write.table(gff,
              file = output,
              quote = F,
              row.names = F,
              col.names = F,
              sep = "\t")
}
