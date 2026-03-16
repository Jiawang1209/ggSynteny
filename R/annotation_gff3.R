#' Read gene and transcript features from a GFF3 file
#'
#' @param gff3 Character path to a GFF3 file.
#' @param gene_types Feature types treated as genes.
#' @param transcript_types Feature types treated as transcripts.
#' @param feature_types Feature types retained for gene structure plotting.
#'
#' @returns A list with `genes`, `transcripts`, and `features` tables.
#' @export
#'
#' @examples NULL
read_gff3_annotation <- function(gff3,
                                 gene_types = c("gene"),
                                 transcript_types = c("mRNA", "transcript"),
                                 feature_types = c("exon", "CDS")) {
  gff_tbl <- readr::read_delim(
    file = gff3,
    delim = "\t",
    comment = "#",
    col_names = FALSE,
    show_col_types = FALSE
  ) |>
    dplyr::filter(nchar(.data$X1) > 0) |>
    dplyr::transmute(
      seqid = as.character(.data$X1),
      source = as.character(.data$X2),
      type = as.character(.data$X3),
      start = as.numeric(.data$X4),
      end = as.numeric(.data$X5),
      score = dplyr::na_if(as.character(.data$X6), "."),
      strand = dplyr::na_if(as.character(.data$X7), "."),
      phase = dplyr::na_if(as.character(.data$X8), "."),
      attributes = as.character(.data$X9)
    ) |>
    dplyr::mutate(
      feature_id = extract_gff3_attribute(.data$attributes, "ID"),
      parent_id = extract_gff3_attribute(.data$attributes, "Parent"),
      gene_name = extract_gff3_attribute(.data$attributes, "Name"),
      transcript_name = extract_gff3_attribute(.data$attributes, "Name"),
      transcript_tag = extract_gff3_attribute(.data$attributes, "tag"),
      canonical_flag = detect_canonical_transcript(.data$attributes)
    )

  genes_tbl <- gff_tbl |>
    dplyr::filter(.data$type %in% gene_types) |>
    dplyr::transmute(
      gene_id = .data$feature_id,
      chr_id = .data$seqid,
      start = .data$start,
      end = .data$end,
      strand = .data$strand,
      gene_name = dplyr::coalesce(.data$gene_name, .data$feature_id)
    )

  transcripts_tbl <- gff_tbl |>
    dplyr::filter(.data$type %in% transcript_types) |>
    dplyr::transmute(
      transcript_id = .data$feature_id,
      gene_id = .data$parent_id,
      chr_id = .data$seqid,
      start = .data$start,
      end = .data$end,
      strand = .data$strand,
      transcript_name = dplyr::coalesce(.data$transcript_name, .data$feature_id),
      transcript_tag = .data$transcript_tag,
      canonical_flag = .data$canonical_flag
    )

  features_tbl <- gff_tbl |>
    dplyr::filter(.data$type %in% feature_types) |>
    dplyr::transmute(
      feature_id = .data$feature_id,
      transcript_id = .data$parent_id,
      chr_id = .data$seqid,
      start = .data$start,
      end = .data$end,
      strand = .data$strand,
      feature_type = .data$type
    ) |>
    dplyr::left_join(
      transcripts_tbl |>
        dplyr::select(transcript_id, gene_id),
      by = "transcript_id"
    )

  cds_lengths <- features_tbl |>
    dplyr::filter(.data$feature_type == "CDS") |>
    dplyr::mutate(width = .data$end - .data$start + 1) |>
    dplyr::group_by(.data$transcript_id) |>
    dplyr::summarise(cds_width = sum(.data$width), .groups = "drop")

  transcripts_tbl <- transcripts_tbl |>
    dplyr::mutate(width = .data$end - .data$start + 1) |>
    dplyr::left_join(cds_lengths, by = "transcript_id") |>
    dplyr::mutate(cds_width = dplyr::coalesce(.data$cds_width, 0))

  list(
    genes = genes_tbl,
    transcripts = transcripts_tbl,
    features = features_tbl
  )
}

#' Attach GFF3 annotation to a ggsynteny data object
#'
#' @param data A `ggsynteny_data` object.
#' @param gff3 Character path to a GFF3 file.
#' @param genome_id Genome ID in `data$genomes` to enrich.
#' @param gene_id_map Optional named character vector mapping GFF3 gene IDs to IDs in `data$genes`.
#' @param transcript_mode One of `"longest_transcript"`, `"longest_cds"`, `"canonical"`,
#'   `"all"`, or `"custom"` for updating gene-level preferred transcripts.
#' @param transcript_ids Optional named character vector mapping `gene_id` to `transcript_id`
#'   when `transcript_mode = "custom"`.
#'
#' @returns An updated `ggsynteny_data` object.
#' @export
#'
#' @examples NULL
add_gff3_annotation <- function(data,
                                gff3,
                                genome_id,
                                gene_id_map = NULL,
                                transcript_mode = "longest_transcript",
                                transcript_ids = NULL) {
  data <- validate_ggsynteny_data(data)
  transcript_mode <- normalize_transcript_mode(transcript_mode)

  if (!genome_id %in% data$genomes$genome_id) {
    stop("`genome_id` must exist in `data$genomes`.", call. = FALSE)
  }

  annotation <- read_gff3_annotation(gff3)
  target_gene_ids <- data$genes |>
    dplyr::filter(.data$genome_id == .env$genome_id) |>
    dplyr::pull(.data$gene_id) |>
    unique()

  gene_lookup <- build_gene_lookup(annotation$genes, gene_id_map) |>
    dplyr::filter(.data$gene_id %in% .env$target_gene_ids)

  transcript_lookup <- build_transcript_lookup(annotation$transcripts, gene_id_map) |>
    dplyr::filter(.data$gene_id %in% .env$target_gene_ids)

  preferred_transcripts <- resolve_preferred_transcripts(
    transcripts = transcript_lookup,
    transcript_mode = transcript_mode,
    transcript_ids = transcript_ids
  )

  genes_annotated <- data$genes |>
    dplyr::left_join(
      gene_lookup |>
        dplyr::select(gene_id, strand_annotation = strand, gene_name_annotation = gene_name),
      by = "gene_id"
    ) |>
    dplyr::left_join(
      preferred_transcripts |>
        dplyr::select(gene_id, transcript_id_annotation = transcript_id),
      by = "gene_id"
    ) |>
    dplyr::mutate(
      strand = dplyr::if_else(
        .data$genome_id == .env$genome_id & !is.na(.data$strand_annotation),
        .data$strand_annotation,
        .data$strand
      ),
      gene_name = dplyr::if_else(
        .data$genome_id == .env$genome_id & !is.na(.data$gene_name_annotation),
        .data$gene_name_annotation,
        dplyr::coalesce(.data$gene_name, .data$gene_id)
      ),
      transcript_id = dplyr::if_else(
        .data$genome_id == .env$genome_id & !is.na(.data$transcript_id_annotation),
        .data$transcript_id_annotation,
        .data$transcript_id
      )
    ) |>
    dplyr::select(-dplyr::any_of(c("strand_annotation", "gene_name_annotation", "transcript_id_annotation")))

  features_tbl <- annotation$features |>
    dplyr::left_join(
      gene_lookup |>
        dplyr::mutate(mapped_gene_id = .data$gene_id) |>
        dplyr::select(gene_id, mapped_gene_id, gene_name),
      by = "gene_id"
    ) |>
    dplyr::transmute(
      genome_id = genome_id,
      gene_id = .data$mapped_gene_id,
      transcript_id = .data$transcript_id,
      chr_id = .data$chr_id,
      start = .data$start,
      end = .data$end,
      strand = .data$strand,
      feature_type = .data$feature_type
    ) |>
    dplyr::filter(!is.na(.data$gene_id), .data$gene_id %in% .env$target_gene_ids)

  data$genes <- genes_annotated
  existing_features <- if ("features" %in% names(data)) data$features else tibble::tibble()
  data$features <- dplyr::bind_rows(existing_features, features_tbl)
  existing_transcripts <- if ("transcripts" %in% names(data)) data$transcripts else tibble::tibble()
  data$transcripts <- dplyr::bind_rows(existing_transcripts, transcript_lookup)

  validate_ggsynteny_data(data)
}

extract_gff3_attribute <- function(attributes, key) {
  pattern <- stringr::str_c("(^|;)", key, "=([^;]+)")
  values <- stringr::str_match(attributes, pattern)[, 3]
  dplyr::na_if(values, "")
}

detect_canonical_transcript <- function(attributes) {
  canonical_patterns <- c(
    "(^|;)canonical=1($|;)",
    "(^|;)canonical=true($|;)",
    "(^|;)is_canonical=true($|;)",
    "MANE_Select",
    "appris_principal",
    "(^|;)tag=canonical($|;)"
  )

  purrr::map_lgl(
    attributes,
    function(x) any(stringr::str_detect(x, canonical_patterns))
  )
}

build_gene_lookup <- function(genes_tbl, gene_id_map = NULL) {
  if (is.null(gene_id_map)) {
    return(genes_tbl)
  }

  tibble::tibble(
    gene_id = unname(gene_id_map),
    gff_gene_id = names(gene_id_map)
  ) |>
    dplyr::left_join(genes_tbl, by = c("gff_gene_id" = "gene_id")) |>
    dplyr::transmute(
      gene_id = .data$gene_id,
      chr_id = .data$chr_id,
      start = .data$start,
      end = .data$end,
      strand = .data$strand,
      gene_name = .data$gene_name
    )
}

build_transcript_lookup <- function(transcripts_tbl, gene_id_map = NULL) {
  if (is.null(gene_id_map)) {
    return(
      transcripts_tbl |>
        dplyr::mutate(width = .data$end - .data$start + 1)
    )
  }

  tibble::tibble(
    gene_id = unname(gene_id_map),
    gff_gene_id = names(gene_id_map)
  ) |>
    dplyr::left_join(transcripts_tbl, by = c("gff_gene_id" = "gene_id")) |>
    dplyr::transmute(
      transcript_id = .data$transcript_id,
      gene_id = .data$gene_id,
      chr_id = .data$chr_id,
      start = .data$start,
      end = .data$end,
      strand = .data$strand,
      transcript_name = .data$transcript_name,
      transcript_tag = .data$transcript_tag,
      canonical_flag = .data$canonical_flag,
      width = .data$width,
      cds_width = .data$cds_width
    )
}

normalize_transcript_mode <- function(transcript_mode) {
  if (length(transcript_mode) != 1) {
    stop("`transcript_mode` must be a single value.", call. = FALSE)
  }

  aliases <- c(
    longest = "longest_transcript",
    cds = "longest_cds"
  )

  if (transcript_mode %in% names(aliases)) {
    transcript_mode <- aliases[[transcript_mode]]
  }

  rlang::arg_match(
    transcript_mode,
    values = c("longest_transcript", "longest_cds", "canonical", "all", "custom")
  )
}

resolve_preferred_transcripts <- function(transcripts,
                                          transcript_mode = "longest_transcript",
                                          transcript_ids = NULL) {
  transcript_mode <- normalize_transcript_mode(transcript_mode)

  select_transcripts(
    transcripts = transcripts,
    transcript_mode = transcript_mode,
    transcript_ids = transcript_ids
  ) |>
    dplyr::group_by(.data$gene_id) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::select(gene_id, transcript_id)
}

select_transcripts <- function(transcripts,
                               transcript_mode = "longest_transcript",
                               transcript_ids = NULL) {
  transcript_mode <- normalize_transcript_mode(transcript_mode)

  if (nrow(transcripts) == 0) {
    return(tibble::tibble(gene_id = character(), transcript_id = character()))
  }

  if (transcript_mode == "all") {
    return(transcripts)
  }

  if (transcript_mode == "custom") {
    if (is.null(transcript_ids) || is.null(names(transcript_ids))) {
      stop("`transcript_ids` must be a named character vector when `transcript_mode = 'custom'`.", call. = FALSE)
    }

    selected <- tibble::tibble(
      gene_id = names(transcript_ids),
      transcript_id = unname(transcript_ids)
    )

    return(
      transcripts |>
        dplyr::inner_join(selected, by = c("gene_id", "transcript_id"))
    )
  }

  if (transcript_mode == "canonical") {
    canonical_hits <- transcripts |>
      dplyr::filter(.data$canonical_flag %in% TRUE)

    if (nrow(canonical_hits) > 0) {
      return(
        canonical_hits |>
          dplyr::arrange(.data$gene_id, dplyr::desc(.data$cds_width), dplyr::desc(.data$width), .data$transcript_id)
      )
    }
  }

  if (transcript_mode == "longest_cds") {
    return(
      transcripts |>
        dplyr::arrange(.data$gene_id, dplyr::desc(.data$cds_width), dplyr::desc(.data$width), .data$transcript_id)
    )
  }

  transcripts |>
    dplyr::arrange(.data$gene_id, dplyr::desc(.data$width), dplyr::desc(.data$cds_width), .data$transcript_id)
}
