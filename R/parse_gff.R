#' Extract chromosome lengths from a GFF3 annotation file
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
parse_gff <- function(gffpath,output){
  
  Chr <- readr::read_delim(file = gffpath, delim = "\t", comment = "#", col_names = F) %>%
    dplyr::filter(X3 == "gene") %>%
    dplyr::select(X1, X4, X5, X9) %>%
    dplyr::mutate(X9 = stringr::str_remove(X9, pattern = ";.*")) %>%
    dplyr::mutate(X9 = stringr::str_remove(X9, pattern = "ID=")) %>%
    dplyr::group_by(X1) %>%
    dplyr::summarise(max(X5)) %>%
    purrr::set_names(c("Chr", "Length"))
  
  utils::write.table(Chr,
                     file = output,
                     quote = FALSE,
                     row.names = FALSE,
                     col.names = FALSE,
                     sep = "\t")
}
