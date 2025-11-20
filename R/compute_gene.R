# 1) 计算线性染色体坐标
compute_chr_layout <- function(chr_df, gap_size = 1e6) {
  chr_df %>%
    mutate(chr_num = stringr::str_extract(Chr, "\\d+") |> as.integer()) %>%
    arrange(chr_num) %>%
    mutate(
      offset       = dplyr::lag(cumsum(Length + gap_size), default = 0),
      global_start = offset,
      global_end   = offset + Length
    )
}

# 2) 把基因坐标映射到全基因组
map_genes_to_linear <- function(genes_df, chr_layout) {
  genes_df %>%
    dplyr::left_join(chr_layout %>% dplyr::select(Chr, offset), by = "Chr") %>%
    dplyr::mutate(
      global_start = offset + Start,
      global_end   = offset + End
    )
}
