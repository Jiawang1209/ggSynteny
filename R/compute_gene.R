#' Compute linear chromosome coordinates
#'
#' @param chr_df A data frame with `genome_id`, `chr_id`, and `length` columns.
#' @param gap_size Numeric gap inserted between chromosomes.
#'
#' @returns A tibble with linear chromosome coordinates.
#' @export
#'
#' @examples NULL
compute_chr_layout <- function(chr_df, gap_size = 1e6) {
  chr_df |>
    dplyr::group_by(.data$genome_id) |>
    dplyr::arrange(order_chromosomes(.data$chr_id), .by_group = TRUE) |>
    dplyr::mutate(
      chr_order = dplyr::row_number(),
      offset = dplyr::lag(cumsum(.data$length + gap_size), default = 0),
      global_start = .data$offset,
      global_end = .data$offset + .data$length
    ) |>
    dplyr::ungroup()
}

#' Map genomic intervals to linear chromosome coordinates
#'
#' @param genes_df A data frame with `genome_id`, `chr_id`, `start`, and `end` columns.
#' @param chr_layout Output of [compute_chr_layout()].
#'
#' @returns A tibble with linear coordinates.
#' @export
#'
#' @examples NULL
map_genes_to_linear <- function(genes_df, chr_layout) {
  genes_df |>
    dplyr::select(-dplyr::any_of(c("offset", "global_start", "global_end"))) |>
    dplyr::left_join(
      chr_layout |> dplyr::select(genome_id, chr_id, offset),
      by = c("genome_id", "chr_id")
    ) |>
    dplyr::mutate(
      global_start = .data$offset + .data$start,
      global_end = .data$offset + .data$end
    )
}

#' Compute linear block coordinates
#'
#' @param block_df A block table produced by [build_blocks_from_anchors()].
#' @param chr_layout Output of [compute_chr_layout()].
#'
#' @returns A tibble with linear coordinates for synteny blocks.
#' @export
#'
#' @examples NULL
compute_block_layout <- function(block_df, chr_layout) {
  chr_offsets <- chr_layout |>
    dplyr::select(genome_id, chr_id, offset)

  block_df |>
    dplyr::left_join(chr_offsets, by = c("genome_a" = "genome_id", "chr_a" = "chr_id")) |>
    dplyr::rename(offset_a = offset) |>
    dplyr::left_join(chr_offsets, by = c("genome_b" = "genome_id", "chr_b" = "chr_id")) |>
    dplyr::rename(offset_b = offset) |>
    dplyr::mutate(
      global_start_a = .data$offset_a + .data$start_a,
      global_end_a = .data$offset_a + .data$end_a,
      global_mid_a = (.data$global_start_a + .data$global_end_a) / 2,
      global_start_b = .data$offset_b + .data$start_b,
      global_end_b = .data$offset_b + .data$end_b,
      global_mid_b = (.data$global_start_b + .data$global_end_b) / 2
    )
}

order_chromosomes <- function(chr) {
  numeric_part <- stringr::str_extract(chr, "[0-9]+")
  has_number <- !is.na(numeric_part)
  dplyr::if_else(
    has_number,
    stringr::str_c("1_", stringr::str_pad(numeric_part, width = 6, pad = "0")),
    stringr::str_c("2_", chr)
  )
}
