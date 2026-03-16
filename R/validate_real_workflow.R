#' Run an end-to-end synteny validation workflow
#'
#' @param mcscanx_gff Character path to an MCScanX gff file.
#' @param collinearity Character path to an MCScanX `.collinearity` file.
#' @param genome1 Label for the first genome.
#' @param genome2 Label for the second genome.
#' @param gff3_genome1 Optional GFF3 annotation for `genome1`.
#' @param gff3_genome2 Optional GFF3 annotation for `genome2`.
#' @param gene_id_map1 Optional named character vector mapping `gff3_genome1` gene IDs to `mcscanx_gff` gene IDs.
#' @param gene_id_map2 Optional named character vector mapping `gff3_genome2` gene IDs to `mcscanx_gff` gene IDs.
#' @param block_id Optional block identifier for micro-synteny plotting. If `NULL`, the first block is used.
#' @param gap_size Numeric gap inserted between chromosomes in genome-level plotting.
#' @param transcript_mode Transcript selection mode passed to annotation and micro-synteny plotting.
#' @param min_pairs Minimum number of anchors required in genome-level plotting.
#' @param flank_genes Number of flanking genes to show in the micro-synteny plot.
#'
#' @returns A list containing standardized data, selected block ID, and ggplot objects.
#' @export
#'
#' @examples NULL
validate_real_synteny_workflow <- function(mcscanx_gff,
                                           collinearity,
                                           genome1 = "Genome1",
                                           genome2 = "Genome2",
                                           gff3_genome1 = NULL,
                                           gff3_genome2 = NULL,
                                           gene_id_map1 = NULL,
                                           gene_id_map2 = NULL,
                                           block_id = NULL,
                                           gap_size = 1e6,
                                           transcript_mode = "longest_transcript",
                                           min_pairs = 5,
                                           flank_genes = 2) {
  transcript_mode <- normalize_transcript_mode(transcript_mode)

  data <- prepare_mcscanx_data(
    synteny = collinearity,
    gff = mcscanx_gff,
    genome1 = genome1,
    genome2 = genome2,
    gap_size = gap_size
  )

  if (!is.null(gff3_genome1)) {
    data <- add_gff3_annotation(
      data = data,
      gff3 = gff3_genome1,
      genome_id = genome1,
      gene_id_map = gene_id_map1,
      transcript_mode = transcript_mode
    )
  }

  if (!is.null(gff3_genome2)) {
    data <- add_gff3_annotation(
      data = data,
      gff3 = gff3_genome2,
      genome_id = genome2,
      gene_id_map = gene_id_map2,
      transcript_mode = transcript_mode
    )
  }

  available_blocks <- data$blocks |>
    dplyr::arrange(dplyr::desc(.data$n_anchors), .data$block_id) |>
    dplyr::pull(.data$block_id)

  if (length(available_blocks) == 0) {
    stop("No synteny blocks were parsed from the input files.", call. = FALSE)
  }

  selected_block_id <- block_id %||% available_blocks[[1]]
  if (!selected_block_id %in% available_blocks) {
    stop("`block_id` must exist in the parsed synteny blocks.", call. = FALSE)
  }

  genome_plot <- plot_Synteny(
    data = data,
    min_pairs = min_pairs
  )

  micro_plot <- plot_microSynteny(
    data = data,
    block_id = selected_block_id,
    flank_genes = flank_genes,
    transcript_mode = transcript_mode,
    show_features = "features" %in% names(data),
    show_gene_label = TRUE
  )

  list(
    data = data,
    block_id = selected_block_id,
    genome_plot = genome_plot,
    micro_plot = micro_plot
  )
}

#' Run an end-to-end jcvi synteny validation workflow
#'
#' @param anchors Character path to a jcvi anchors file.
#' @param bed_a Character path to the first jcvi BED file.
#' @param bed_b Character path to the second jcvi BED file.
#' @param genome1 Label for the first genome.
#' @param genome2 Label for the second genome.
#' @param gff3_genome1 Optional GFF3 annotation for `genome1`.
#' @param gff3_genome2 Optional GFF3 annotation for `genome2`.
#' @param gene_id_map1 Optional named character vector mapping `gff3_genome1` gene IDs to jcvi gene IDs.
#' @param gene_id_map2 Optional named character vector mapping `gff3_genome2` gene IDs to jcvi gene IDs.
#' @param block_id Optional block identifier for micro-synteny plotting. If `NULL`, the first block is used.
#' @param gap_size Numeric gap inserted between chromosomes in genome-level plotting.
#' @param transcript_mode Transcript selection mode passed to annotation and micro-synteny plotting.
#' @param min_pairs Minimum number of anchors required in genome-level plotting.
#' @param flank_genes Number of flanking genes to show in the micro-synteny plot.
#'
#' @returns A list containing standardized data, selected block ID, and ggplot objects.
#' @export
#'
#' @examples NULL
validate_jcvi_synteny_workflow <- function(anchors,
                                           bed_a,
                                           bed_b,
                                           genome1 = "Genome1",
                                           genome2 = "Genome2",
                                           gff3_genome1 = NULL,
                                           gff3_genome2 = NULL,
                                           gene_id_map1 = NULL,
                                           gene_id_map2 = NULL,
                                           block_id = NULL,
                                           gap_size = 1e6,
                                           transcript_mode = "longest_transcript",
                                           min_pairs = 5,
                                           flank_genes = 2) {
  transcript_mode <- normalize_transcript_mode(transcript_mode)

  data <- prepare_jcvi_data(
    anchors = anchors,
    bed_a = bed_a,
    bed_b = bed_b,
    genome1 = genome1,
    genome2 = genome2,
    gap_size = gap_size
  )

  if (!is.null(gff3_genome1)) {
    data <- add_gff3_annotation(
      data = data,
      gff3 = gff3_genome1,
      genome_id = genome1,
      gene_id_map = gene_id_map1,
      transcript_mode = transcript_mode
    )
  }

  if (!is.null(gff3_genome2)) {
    data <- add_gff3_annotation(
      data = data,
      gff3 = gff3_genome2,
      genome_id = genome2,
      gene_id_map = gene_id_map2,
      transcript_mode = transcript_mode
    )
  }

  available_blocks <- data$blocks |>
    dplyr::arrange(dplyr::desc(.data$n_anchors), .data$block_id) |>
    dplyr::pull(.data$block_id)

  if (length(available_blocks) == 0) {
    stop("No synteny blocks were parsed from the input files.", call. = FALSE)
  }

  selected_block_id <- block_id %||% available_blocks[[1]]
  if (!selected_block_id %in% available_blocks) {
    stop("`block_id` must exist in the parsed synteny blocks.", call. = FALSE)
  }

  genome_plot <- plot_Synteny(
    data = data,
    min_pairs = min_pairs
  )

  micro_plot <- plot_microSynteny(
    data = data,
    block_id = selected_block_id,
    flank_genes = flank_genes,
    transcript_mode = transcript_mode,
    show_features = "features" %in% names(data),
    show_gene_label = TRUE
  )

  list(
    data = data,
    block_id = selected_block_id,
    genome_plot = genome_plot,
    micro_plot = micro_plot
  )
}
