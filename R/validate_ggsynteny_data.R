#' Validate a ggsynteny data object
#'
#' @param data A `ggsynteny_data` object.
#'
#' @returns The validated input object.
#' @export
#'
#' @examples NULL
validate_ggsynteny_data <- function(data) {
  if (!is.list(data)) {
    stop("`data` must be a list-like object.", call. = FALSE)
  }

  required_tables <- c("genomes", "chromosomes", "genes", "anchors", "blocks")
  missing_tables <- setdiff(required_tables, names(data))
  if (length(missing_tables) > 0) {
    stop(
      "Missing required tables in `ggsynteny_data`: ",
      paste(missing_tables, collapse = ", "),
      call. = FALSE
    )
  }

  check_required_columns(data$genomes, c("genome_id", "genome_label", "order"), "genomes")
  check_required_columns(
    data$chromosomes,
    c("genome_id", "chr_id", "length", "chr_order", "global_start", "global_end"),
    "chromosomes"
  )
  check_required_columns(
    data$genes,
    c("genome_id", "gene_id", "chr_id", "start", "end"),
    "genes"
  )
  check_required_columns(
    data$anchors,
    c("anchor_id", "block_id", "genome_a", "gene_a", "chr_a", "genome_b", "gene_b", "chr_b"),
    "anchors"
  )
  check_required_columns(
    data$blocks,
    c("block_id", "genome_a", "chr_a", "start_a", "end_a", "genome_b", "chr_b", "start_b", "end_b"),
    "blocks"
  )

  known_genomes <- unique(data$genomes$genome_id)
  referenced_genomes <- unique(c(
    data$chromosomes$genome_id,
    data$genes$genome_id,
    data$anchors$genome_a,
    data$anchors$genome_b,
    data$blocks$genome_a,
    data$blocks$genome_b
  ))

  if (!all(referenced_genomes %in% known_genomes)) {
    stop("All tables must reference genome IDs defined in `genomes`.", call. = FALSE)
  }

  chromosome_keys <- paste(data$chromosomes$genome_id, data$chromosomes$chr_id, sep = "::")
  gene_chromosome_keys <- paste(data$genes$genome_id, data$genes$chr_id, sep = "::")
  if (!all(gene_chromosome_keys %in% chromosome_keys)) {
    stop("Each gene must map to a chromosome listed in `chromosomes`.", call. = FALSE)
  }

  block_ids <- unique(data$blocks$block_id)
  if (!all(unique(data$anchors$block_id) %in% block_ids)) {
    stop("Each anchor must map to a block listed in `blocks`.", call. = FALSE)
  }

  if ("features" %in% names(data) && nrow(data$features) > 0) {
    check_required_columns(
      data$features,
      c("genome_id", "gene_id", "chr_id", "start", "end", "feature_type"),
      "features"
    )

    gene_keys <- paste(data$genes$genome_id, data$genes$gene_id, sep = "::")
    feature_gene_keys <- paste(data$features$genome_id, data$features$gene_id, sep = "::")
    if (!all(feature_gene_keys %in% gene_keys)) {
      stop("Each feature must map to a gene listed in `genes`.", call. = FALSE)
    }
  }

  if ("transcripts" %in% names(data) && nrow(data$transcripts) > 0) {
    check_required_columns(
      data$transcripts,
      c("gene_id", "transcript_id", "chr_id", "start", "end"),
      "transcripts"
    )
  }

  if (is.null(class(data)) || !"ggsynteny_data" %in% class(data)) {
    class(data) <- unique(c("ggsynteny_data", class(data)))
  }

  data
}

check_required_columns <- function(tbl, required_columns, table_name) {
  missing_columns <- setdiff(required_columns, names(tbl))
  if (length(missing_columns) > 0) {
    stop(
      "Table `", table_name, "` is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}
