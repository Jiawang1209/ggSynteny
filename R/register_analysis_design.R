#' Register genomes for a multi-genome synteny analysis
#'
#' @param genome_id Character vector of genome identifiers.
#' @param annotation Optional character vector of annotation file paths.
#' @param protein Optional character vector of protein FASTA paths.
#' @param cds Optional character vector of CDS FASTA paths.
#' @param genome_label Optional display labels. Defaults to `genome_id`.
#' @param order Optional integer order. Defaults to input order.
#' @param group Optional character vector for grouping genomes.
#'
#' @returns A tibble with one row per genome.
#' @export
#'
#' @examples NULL
register_genomes <- function(genome_id,
                             annotation = NULL,
                             protein = NULL,
                             cds = NULL,
                             genome_label = NULL,
                             order = NULL,
                             group = NULL) {
  genome_id <- as.character(genome_id)

  if (length(genome_id) == 0) {
    stop("`genome_id` must contain at least one genome.", call. = FALSE)
  }

  if (anyNA(genome_id) || any(genome_id == "")) {
    stop("`genome_id` must not contain missing or empty values.", call. = FALSE)
  }

  if (anyDuplicated(genome_id)) {
    stop("`genome_id` values must be unique.", call. = FALSE)
  }

  n <- length(genome_id)

  annotation <- normalize_registry_field(annotation, n, "annotation")
  protein <- normalize_registry_field(protein, n, "protein")
  cds <- normalize_registry_field(cds, n, "cds")
  genome_label <- normalize_registry_field(genome_label, n, "genome_label", default = genome_id)
  group <- normalize_registry_field(group, n, "group")

  if (is.null(order)) {
    order <- seq_len(n)
  }

  if (length(order) != n || anyNA(order) || anyDuplicated(order)) {
    stop("`order` must have one unique non-missing value per genome.", call. = FALSE)
  }

  tibble::tibble(
    genome_id = genome_id,
    genome_label = genome_label,
    annotation_path = annotation,
    protein_fasta = protein,
    cds_fasta = cds,
    group = group,
    order = as.integer(order)
  ) |>
    dplyr::arrange(.data$order)
}

#' Define ordered pairwise comparisons for a multi-genome analysis
#'
#' @param genomes A genome registry produced by [register_genomes()].
#' @param query Character vector of query genome IDs.
#' @param subject Character vector of subject genome IDs.
#' @param comparison_id Optional character vector of comparison IDs.
#' @param order Optional integer order. Defaults to input order.
#' @param analysis_tool Character vector identifying the downstream tool.
#'
#' @returns A tibble with one row per comparison.
#' @export
#'
#' @examples NULL
define_comparisons <- function(genomes,
                               query,
                               subject,
                               comparison_id = NULL,
                               order = NULL,
                               analysis_tool = "mcscanx") {
  genomes <- validate_genome_registry(genomes)
  query <- as.character(query)
  subject <- as.character(subject)

  if (length(query) == 0 || length(subject) == 0 || length(query) != length(subject)) {
    stop("`query` and `subject` must have the same non-zero length.", call. = FALSE)
  }

  known_genomes <- genomes$genome_id
  if (!all(query %in% known_genomes)) {
    stop("All `query` genomes must be registered in `genomes`.", call. = FALSE)
  }
  if (!all(subject %in% known_genomes)) {
    stop("All `subject` genomes must be registered in `genomes`.", call. = FALSE)
  }

  n <- length(query)
  analysis_tool <- normalize_registry_field(analysis_tool, n, "analysis_tool", default = "mcscanx")

  if (is.null(comparison_id)) {
    comparison_id <- stringr::str_c(query, "vs", subject, sep = "_")
  }

  if (length(comparison_id) != n || anyNA(comparison_id) || any(comparison_id == "")) {
    stop("`comparison_id` must have one non-missing value per comparison.", call. = FALSE)
  }

  if (anyDuplicated(comparison_id)) {
    stop("`comparison_id` values must be unique.", call. = FALSE)
  }

  if (is.null(order)) {
    order <- seq_len(n)
  }

  if (length(order) != n || anyNA(order) || anyDuplicated(order)) {
    stop("`order` must have one unique non-missing value per comparison.", call. = FALSE)
  }

  tibble::tibble(
    comparison_id = comparison_id,
    query_genome = query,
    subject_genome = subject,
    analysis_tool = analysis_tool,
    order = as.integer(order)
  ) |>
    dplyr::arrange(.data$order)
}

#' Create a multi-genome analysis design
#'
#' @param genomes A genome registry produced by [register_genomes()].
#' @param comparisons A comparison registry produced by [define_comparisons()].
#'
#' @returns A list with class `ggsynteny_design`.
#' @export
#'
#' @examples NULL
create_analysis_design <- function(genomes, comparisons) {
  genomes <- validate_genome_registry(genomes)
  comparisons <- validate_comparison_registry(comparisons, genomes = genomes)

  structure(
    list(
      genomes = genomes,
      comparisons = comparisons
    ),
    class = "ggsynteny_design"
  )
}

#' Validate a multi-genome analysis design
#'
#' @param design A design created by [create_analysis_design()].
#'
#' @returns The validated design.
#' @export
#'
#' @examples NULL
validate_analysis_design <- function(design) {
  if (!is.list(design)) {
    stop("`design` must be a list-like object.", call. = FALSE)
  }

  if (!all(c("genomes", "comparisons") %in% names(design))) {
    stop("`design` must contain `genomes` and `comparisons`.", call. = FALSE)
  }

  design$genomes <- validate_genome_registry(design$genomes)
  design$comparisons <- validate_comparison_registry(design$comparisons, genomes = design$genomes)

  if (is.null(class(design)) || !"ggsynteny_design" %in% class(design)) {
    class(design) <- unique(c("ggsynteny_design", class(design)))
  }

  design
}

validate_genome_registry <- function(genomes) {
  required_columns <- c(
    "genome_id", "genome_label", "annotation_path",
    "protein_fasta", "cds_fasta", "group", "order"
  )
  check_required_columns(genomes, required_columns, "genomes")

  if (anyDuplicated(genomes$genome_id)) {
    stop("`genomes$genome_id` must be unique.", call. = FALSE)
  }

  if (anyDuplicated(genomes$order)) {
    stop("`genomes$order` must be unique.", call. = FALSE)
  }

  genomes |>
    dplyr::arrange(.data$order)
}

validate_comparison_registry <- function(comparisons, genomes) {
  required_columns <- c(
    "comparison_id", "query_genome", "subject_genome", "analysis_tool", "order"
  )
  check_required_columns(comparisons, required_columns, "comparisons")

  if (anyDuplicated(comparisons$comparison_id)) {
    stop("`comparisons$comparison_id` must be unique.", call. = FALSE)
  }

  if (anyDuplicated(comparisons$order)) {
    stop("`comparisons$order` must be unique.", call. = FALSE)
  }

  known_genomes <- genomes$genome_id
  if (!all(comparisons$query_genome %in% known_genomes)) {
    stop("All `query_genome` values must exist in `genomes`.", call. = FALSE)
  }
  if (!all(comparisons$subject_genome %in% known_genomes)) {
    stop("All `subject_genome` values must exist in `genomes`.", call. = FALSE)
  }

  comparisons |>
    dplyr::arrange(.data$order)
}

normalize_registry_field <- function(x, n, field_name, default = NULL) {
  if (is.null(x)) {
    return(default %||% rep(NA_character_, n))
  }

  if (length(x) == 1 && n > 1) {
    return(rep(x, n))
  }

  if (length(x) != n) {
    stop("`", field_name, "` must have length ", n, ".", call. = FALSE)
  }

  x
}
