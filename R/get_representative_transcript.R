#' Extract representative protein transcript from a FASTA file
#'
#' @param fasta Character
#' Path to the input amino acid FASTA file to be processed.
#' @param sep Character
#' Separator used in the sequence IDs to distinguish gene ID and transcript ID. 
#' @param output Character
#'  Path to the output FASTA file that will store the representative (longest) transcript per gene.
#'
#' @returns None
#' @export
#'
#' @examples NULL
#' 
get_representative_transcript <- function(fasta, sep, output){

  # load fasta 
  df <- Biostrings::readAAStringSet(filepath = fasta) %>%
    as.data.frame() 
  
  df2 <- df %>%
    tibble::rownames_to_column(var = "ID") %>%
    purrr::set_names(c("ID", "fasta")) %>%
    dplyr::mutate(length = stringr::str_length(fasta)) %>%
    dplyr::mutate(ID2 = stringr::str_remove(ID, pattern = stringr::str_c("\\", sep, ".*", sep = ""))) %>%
    dplyr::arrange(ID2, desc(length)) %>%
    dplyr::distinct(ID2, .keep_all = T) %>%
    dplyr::mutate(out = str_c(">",ID2,"\n",fasta))
  
  df3 <- df2 %>%
    dplyr::select(ID, ID2, fasta, length)
  
  df2 %>%
    dplyr::select(out) %>%
    write.table(file = output,
                quote = F,
                row.names = F,
                col.names = F)

  
  cat("The", fasta, "have", dim(df)[1], "fasta!\n")
  cat("Clean done! \n")
  cat("The clean fasta contains",  dim(df2)[1], "fasta!")
  
}
