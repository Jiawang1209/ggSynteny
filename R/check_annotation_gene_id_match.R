#' Check gene ID compatibility between synteny data and a GFF3 file
#'
#' @param data A `ggsynteny_data` object.
#' @param gff3 Character path to a GFF3 file.
#' @param genome_id Genome ID in `data$genomes` to inspect.
#' @param gene_id_map Optional named character vector mapping GFF3 gene IDs to
#'   gene IDs used in the synteny object.
#'
#' @returns A list summarizing gene ID overlap and unmatched IDs.
#' @export
#'
#' @examples NULL
check_annotation_gene_id_match <- function(data,
                                           gff3,
                                           genome_id,
                                           gene_id_map = NULL) {
  data <- validate_ggsynteny_data(data)

  if (!genome_id %in% data$genomes$genome_id) {
    stop("`genome_id` must exist in `data$genomes`.", call. = FALSE)
  }

  annotation <- read_gff3_annotation(gff3)
  data_gene_ids <- data$genes |>
    dplyr::filter(.data$genome_id == .env$genome_id) |>
    dplyr::pull(.data$gene_id) |>
    unique()

  annotation_gene_ids <- annotation$genes |>
    dplyr::pull(.data$gene_id) |>
    unique()

  mapped_gene_ids <- if (is.null(gene_id_map)) {
    annotation_gene_ids
  } else {
    unname(gene_id_map[names(gene_id_map) %in% annotation_gene_ids])
  }

  matched_ids <- intersect(data_gene_ids, mapped_gene_ids)
  unmatched_data_ids <- setdiff(data_gene_ids, mapped_gene_ids)
  unmatched_annotation_ids <- if (is.null(gene_id_map)) {
    setdiff(annotation_gene_ids, data_gene_ids)
  } else {
    annotation_gene_ids[!annotation_gene_ids %in% names(gene_id_map)[unname(gene_id_map) %in% data_gene_ids]]
  }

  list(
    genome_id = genome_id,
    n_data_genes = length(data_gene_ids),
    n_annotation_genes = length(annotation_gene_ids),
    n_matched_genes = length(matched_ids),
    match_rate = if (length(data_gene_ids) == 0) 0 else length(matched_ids) / length(data_gene_ids),
    matched_gene_ids = matched_ids,
    unmatched_data_gene_ids = unmatched_data_ids,
    unmatched_annotation_gene_ids = unmatched_annotation_ids
  )
}
