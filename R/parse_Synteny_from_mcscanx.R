parse_Synteny_from_mcscanx <- function(synteny, gff){
  
  # synteny = "Genome_data_2/Athaliana_Osativa.collinearity"
  # gff = "Genome_data_2/Athaliana_Osativa.gff"
  
  # load gff file
  gff <- readr::read_delim(gff,
                           delim = "\t",
                           col_names = F)
  
  gff_Chr <- gff %>%
    dplyr::mutate(tmp = ifelse(X3 < X4, X4, X3)) %>%
    dplyr::group_by(X1) %>%
    dplyr::summarise(Length = max(tmp))
  
  # load synteny output
  synteny_df <- readr::read_delim(synteny, 
                                  delim = "\t",
                                  comment = "#",
                                  col_names = F) %>%
    dplyr::mutate(X1 = stringr::str_remove(X1, pattern = "^\\s+")) %>%
    tidyr::separate(col = X1, sep = "-", into = c("block", "ID"), remove = F) %>%
    dplyr::mutate(ID = stringr::str_remove(ID, pattern = "^\\s+")) %>%
    dplyr::mutate(ID = stringr::str_remove(ID, pattern = ":")) %>%
    dplyr::mutate(ID2 = as.numeric(ID))
  
  # get synteny block and location
  synteny_block <- synteny_df %>%
    dplyr::group_by(block) %>%
    dplyr::summarise(min = min(ID2),
                     max = max(ID2)) %>%
    tidyr::pivot_longer(cols = -block,
                        names_to = "Kind",
                        values_to = "ID2")
  
  synteny_block_species1 <- synteny_block %>%
    dplyr::left_join(synteny_df %>% dplyr::select(2,3,4,7), by = c("block", "ID2")) %>%
    dplyr::left_join(gff, by = c("X2")) %>%
    purrr::set_names(c("block","Kind", "ID2", "ID", "Gene", "Chr", "Start", "End")) %>%
    dplyr::group_by(block, Chr) %>%
    dplyr::summarise(min = min(Start),
                     max = max(End)) %>%
    dplyr::ungroup()
  
  synteny_block_species2 <- synteny_block %>%
    dplyr::left_join(synteny_df %>% dplyr::select(2,3,5,7), by = c("block", "ID2")) %>%
    dplyr::left_join(gff, by = c( "X3" = "X2")) %>%
    purrr::set_names(c("block","Kind", "ID2", "ID", "Gene", "Chr", "Start", "End")) %>%
    dplyr::group_by(block, Chr) %>%
    dplyr::summarise(min = min(Start),
                     max = max(End)) %>%
    dplyr::ungroup()
  
  synteny_block_out <- dplyr::left_join(synteny_block_species1, 
                                        synteny_block_species2,
                                        by = c("block")) %>%
    purrr::set_names(c("block", 
                       "Chr1", "Start1", "End1", 
                       "Chr2", "Start2", "End2"))
  
  return(list(synteny_block_out, 
              gff_Chr))
}
