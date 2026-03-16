#' Read a jcvi BED file
#'
#' @param path Character path to a jcvi BED file.
#' @param genome Optional genome label.
#'
#' @returns A tibble with gene coordinates.
#' @export
#'
#' @examples NULL
read_jcvi_bed <- function(path, genome = NA_character_) {
  bed <- readr::read_delim(
    file = path,
    delim = "\t",
    comment = "#",
    col_names = FALSE,
    show_col_types = FALSE,
    trim_ws = TRUE
  )

  if (ncol(bed) < 4) {
    stop("jcvi BED files must have at least 4 columns.", call. = FALSE)
  }

  strand_col <- if (ncol(bed) >= 6) "X6" else NULL

  bed |>
    dplyr::transmute(
      genome_id = genome,
      chr_id = as.character(.data$X1),
      chr_label = as.character(.data$X1),
      gene_id = as.character(.data$X4),
      gene_name = as.character(.data$X4),
      transcript_id = NA_character_,
      start = pmin(as.numeric(.data$X2), as.numeric(.data$X3), na.rm = TRUE),
      end = pmax(as.numeric(.data$X2), as.numeric(.data$X3), na.rm = TRUE),
      strand = if (is.null(strand_col)) NA_character_ else as.character(.data[[strand_col]])
    ) |>
    dplyr::group_by(.data$chr_id) |>
    dplyr::arrange(.data$start, .data$end, .by_group = TRUE) |>
    dplyr::mutate(gene_order = dplyr::row_number()) |>
    dplyr::ungroup()
}

#' Parse synteny blocks from jcvi anchors and BED files
#'
#' @param anchors Character path to a jcvi anchors file.
#' @param bed_a Character path to the first jcvi BED file.
#' @param bed_b Character path to the second jcvi BED file.
#' @param genome1 Label for the upper genome in the plot.
#' @param genome2 Label for the lower genome in the plot.
#' @param gap_size Numeric gap inserted between chromosomes when computing linear coordinates.
#'
#' @returns A standardized `ggsynteny_data` object.
#' @export
#'
#' @examples NULL
prepare_jcvi_data <- function(anchors,
                              bed_a,
                              bed_b,
                              genome1 = "Genome1",
                              genome2 = "Genome2",
                              gap_size = 1e6) {
  genes_a <- read_jcvi_bed(bed_a, genome = genome1)
  genes_b <- read_jcvi_bed(bed_b, genome = genome2)
  anchors_raw <- read_jcvi_anchors(anchors)

  genome_tbl <- tibble::tibble(
    genome_id = c(genome1, genome2),
    genome_label = c(genome1, genome2),
    group = NA_character_,
    order = c(1L, 2L)
  )

  genes_tbl <- dplyr::bind_rows(genes_a, genes_b)

  chromosomes_tbl <- build_mcscanx_chr_table(genes_tbl, genome_tbl) |>
    compute_chr_layout(gap_size = gap_size)

  genes_tbl <- map_genes_to_linear(genes_tbl, chromosomes_tbl)

  anchors_tbl <- build_jcvi_anchor_table(
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
      metadata = list(source = "jcvi", gap_size = gap_size),
      synteny_blocks = blocks_tbl
    ),
    class = "ggsynteny_data"
  )

  validate_ggsynteny_data(data)
}

#' Parse synteny blocks from jcvi output
#'
#' @param anchors Character path to a jcvi anchors file.
#' @param bed_a Character path to the first jcvi BED file.
#' @param bed_b Character path to the second jcvi BED file.
#' @param genome1 Label for the first genome.
#' @param genome2 Label for the second genome.
#' @param gap_size Numeric gap inserted between chromosomes when computing linear coordinates.
#'
#' @returns A legacy list containing block coordinates and chromosome lengths.
#' @export
#'
#' @examples NULL
parse_Synteny_from_jcvi <- function(anchors,
                                    bed_a,
                                    bed_b,
                                    genome1 = "Genome1",
                                    genome2 = "Genome2",
                                    gap_size = 1e6) {
  parsed <- prepare_jcvi_data(
    anchors = anchors,
    bed_a = bed_a,
    bed_b = bed_b,
    genome1 = genome1,
    genome2 = genome2,
    gap_size = gap_size
  )

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

read_jcvi_anchors <- function(path) {
  lines <- readr::read_lines(path)

  if (length(lines) == 0) {
    stop("jcvi anchors file is empty.", call. = FALSE)
  }

  block_index <- -1L
  anchor_rows <- vector("list", length(lines))
  anchor_count <- 0L

  for (line in lines) {
    line <- stringr::str_trim(line)

    if (identical(line, "")) {
      next
    }

    if (stringr::str_starts(line, "###")) {
      block_index <- block_index + 1L
      next
    }

    if (stringr::str_starts(line, "#")) {
      next
    }

    fields <- stringr::str_split(line, "\\s+", simplify = TRUE)
    if (ncol(fields) < 2) {
      next
    }

    if (block_index < 0L) {
      block_index <- 0L
    }

    anchor_count <- anchor_count + 1L
    score_value <- suppressWarnings(as.numeric(fields[1, min(3L, ncol(fields))]))

    anchor_rows[[anchor_count]] <- tibble::tibble(
      block_id = as.character(block_index),
      anchor_rank = anchor_count,
      gene_a = fields[1, 1],
      gene_b = fields[1, 2],
      score = score_value
    )
  }

  anchors_tbl <- dplyr::bind_rows(anchor_rows[seq_len(anchor_count)])

  if (nrow(anchors_tbl) == 0) {
    stop("No anchor pairs were parsed from the jcvi anchors file.", call. = FALSE)
  }

  anchors_tbl |>
    dplyr::group_by(.data$block_id) |>
    dplyr::mutate(anchor_rank = dplyr::row_number()) |>
    dplyr::ungroup()
}

build_jcvi_anchor_table <- function(anchors_raw, genes_tbl, genome_a, genome_b) {
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
      source = "jcvi",
      genome_a = genome_a,
      genome_b = genome_b,
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
