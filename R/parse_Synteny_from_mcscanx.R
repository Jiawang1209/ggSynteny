#' Read an MCScanX gff file
#'
#' @param path Character path to the MCScanX gff file.
#' @param genome Optional genome label.
#'
#' @returns A tibble with gene coordinates.
#' @export
#'
#' @examples NULL
read_mcscanx_gff <- function(path, genome = NA_character_) {
  gff <- readr::read_delim(
    file = path,
    delim = "\t",
    col_names = FALSE,
    show_col_types = FALSE
  )

  if (ncol(gff) < 4) {
    stop("MCScanX gff file must have at least 4 columns.", call. = FALSE)
  }

  gff |>
    dplyr::transmute(
      genome_id = genome,
      chr_id = as.character(.data$X1),
      chr_label = as.character(.data$X1),
      gene_id = as.character(.data$X2),
      gene_name = as.character(.data$X2),
      transcript_id = NA_character_,
      start = pmin(as.numeric(.data$X3), as.numeric(.data$X4), na.rm = TRUE),
      end = pmax(as.numeric(.data$X3), as.numeric(.data$X4), na.rm = TRUE),
      strand = NA_character_
    ) |>
    dplyr::arrange(.data$chr_id, .data$start, .data$end, .data$gene_id)
}

#' Parse synteny blocks from an MCScanX collinearity file
#'
#' @param synteny Character path to the MCScanX `.collinearity` file.
#' @param gff Character path to the MCScanX gff file.
#' @param genome1 Label for the upper genome in the plot.
#' @param genome2 Label for the lower genome in the plot.
#' @param gap_size Numeric gap inserted between chromosomes when computing linear coordinates.
#'
#' @returns A standardized `ggsynteny_data` object.
#' @export
#'
#' @examples NULL
prepare_mcscanx_data <- function(synteny,
                                 gff,
                                 genome1 = "Genome1",
                                 genome2 = "Genome2",
                                 gap_size = 1e6) {
  genes_raw <- read_mcscanx_gff(gff)

  anchors_raw <- readr::read_delim(
    file = synteny,
    delim = "\t",
    comment = "#",
    col_names = FALSE,
    show_col_types = FALSE
  )

  if (ncol(anchors_raw) < 3) {
    stop("MCScanX collinearity file must have at least 3 columns.", call. = FALSE)
  }

  anchors_raw <- anchors_raw |>
    dplyr::mutate(
      X1 = stringr::str_squish(as.character(.data$X1)),
      block_id = stringr::str_extract(.data$X1, "^[0-9]+"),
      anchor_rank = readr::parse_number(.data$X1),
      gene_a = as.character(.data$X2),
      gene_b = as.character(.data$X3)
    ) |>
    dplyr::filter(!is.na(.data$block_id), !is.na(.data$gene_a), !is.na(.data$gene_b)) |>
    dplyr::select(block_id, anchor_rank, gene_a, gene_b)

  genome_tbl <- tibble::tibble(
    genome_id = c(genome1, genome2),
    genome_label = c(genome1, genome2),
    group = NA_character_,
    order = c(1L, 2L)
  )

  genes_a <- genes_raw |>
    dplyr::filter(.data$gene_id %in% anchors_raw$gene_a) |>
    dplyr::mutate(genome_id = genome1) |>
    dplyr::group_by(.data$chr_id) |>
    dplyr::arrange(.data$start, .by_group = TRUE) |>
    dplyr::mutate(gene_order = dplyr::row_number()) |>
    dplyr::ungroup()

  genes_b <- genes_raw |>
    dplyr::filter(.data$gene_id %in% anchors_raw$gene_b) |>
    dplyr::mutate(genome_id = genome2) |>
    dplyr::group_by(.data$chr_id) |>
    dplyr::arrange(.data$start, .by_group = TRUE) |>
    dplyr::mutate(gene_order = dplyr::row_number()) |>
    dplyr::ungroup()

  genes_tbl <- dplyr::bind_rows(genes_a, genes_b)

  chromosomes_tbl <- build_mcscanx_chr_table(genes_tbl, genome_tbl) |>
    compute_chr_layout(gap_size = gap_size)

  genes_tbl <- map_genes_to_linear(genes_tbl, chromosomes_tbl)

  anchors_tbl <- build_mcscanx_anchor_table(
    anchors_raw = anchors_raw,
    genes_tbl = genes_tbl,
    genome_a = genome1,
    genome_b = genome2
  )

  blocks_tbl <- build_blocks_from_anchors(anchors_tbl) |>
    compute_block_layout(chromosomes_tbl)

  data <- structure(
    list(
      genomes = genome_tbl,
      chromosomes = chromosomes_tbl,
      genes = genes_tbl,
      anchors = anchors_tbl,
      blocks = blocks_tbl,
      metadata = list(source = "mcscanx", gap_size = gap_size),
      synteny_blocks = blocks_tbl
    ),
    class = "ggsynteny_data"
  )

  validate_ggsynteny_data(data)
}

#' Parse synteny blocks from MCScanX output
#'
#' @param synteny Character path to the MCScanX `.collinearity` file.
#' @param gff Character path to the MCScanX gff file.
#'
#' @returns A legacy list containing block coordinates and chromosome lengths.
#' @export
#'
#' @examples NULL
parse_Synteny_from_mcscanx <- function(synteny, gff) {
  parsed <- prepare_mcscanx_data(synteny = synteny, gff = gff)

  synteny_block_out <- parsed$blocks |>
    dplyr::transmute(
      block = .data$block_id,
      Chr1 = .data$chr_a,
      Start1 = .data$start_a,
      End1 = .data$end_a,
      Chr2 = .data$chr_b,
      Start2 = .data$start_b,
      End2 = .data$end_b
    )

  gff_Chr <- parsed$chromosomes |>
    dplyr::distinct(chr_id, length) |>
    dplyr::transmute(Chr = .data$chr_id, Length = .data$length)

  list(synteny_block_out, gff_Chr)
}

build_mcscanx_chr_table <- function(genes_tbl, genome_tbl) {
  genes_tbl |>
    dplyr::group_by(.data$genome_id, .data$chr_id, .data$chr_label) |>
    dplyr::summarise(length = max(.data$end, na.rm = TRUE), .groups = "drop") |>
    dplyr::left_join(genome_tbl |> dplyr::select(genome_id, order), by = "genome_id") |>
    dplyr::mutate(panel = dplyr::if_else(.data$order == 1L, "top", "bottom"))
}

build_mcscanx_anchor_table <- function(anchors_raw, genes_tbl, genome_a, genome_b) {
  genes_a <- genes_tbl |>
    dplyr::filter(.data$genome_id == genome_a) |>
    dplyr::select(
      gene_a = gene_id,
      chr_a = chr_id,
      start_a = start,
      end_a = end,
      global_start_a = global_start,
      global_end_a = global_end
    )

  genes_b <- genes_tbl |>
    dplyr::filter(.data$genome_id == genome_b) |>
    dplyr::select(
      gene_b = gene_id,
      chr_b = chr_id,
      start_b = start,
      end_b = end,
      global_start_b = global_start,
      global_end_b = global_end
    )

  anchors_raw |>
    dplyr::left_join(genes_a, by = "gene_a") |>
    dplyr::left_join(genes_b, by = "gene_b") |>
    dplyr::filter(
      !is.na(.data$chr_a),
      !is.na(.data$chr_b),
      !is.na(.data$start_a),
      !is.na(.data$start_b)
    ) |>
    dplyr::mutate(
      anchor_id = dplyr::row_number(),
      source = "mcscanx",
      genome_a = genome_a,
      genome_b = genome_b,
      score = NA_real_,
      evalue = NA_real_
    ) |>
    dplyr::select(
      anchor_id,
      block_id,
      source,
      anchor_rank,
      genome_a,
      chr_a,
      gene_a,
      start_a,
      end_a,
      global_start_a,
      global_end_a,
      genome_b,
      chr_b,
      gene_b,
      start_b,
      end_b,
      global_start_b,
      global_end_b,
      score,
      evalue
    )
}

#' Build synteny blocks from anchor pairs
#'
#' @param anchors A standardized anchor table.
#'
#' @returns A standardized block table.
#' @export
#'
#' @examples NULL
build_blocks_from_anchors <- function(anchors) {
  anchors |>
    dplyr::group_by(.data$block_id, .data$source, .data$genome_a, .data$genome_b) |>
    dplyr::summarise(
      chr_a = dplyr::first(.data$chr_a),
      start_a = min(.data$start_a, na.rm = TRUE),
      end_a = max(.data$end_a, na.rm = TRUE),
      chr_b = dplyr::first(.data$chr_b),
      start_b = min(.data$start_b, na.rm = TRUE),
      end_b = max(.data$end_b, na.rm = TRUE),
      n_anchors = dplyr::n(),
      block_score = NA_real_,
      .groups = "drop"
    )
}
