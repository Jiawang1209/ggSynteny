#' Prepare ordered MCScanX inputs from a multi-genome analysis design
#'
#' @param design A design created by [create_analysis_design()].
#' @param outdir Character output directory for prepared files.
#' @param protein_sep Optional separator used to collapse transcript-level protein
#'   FASTA IDs into gene-level representative sequences.
#' @param prepare_protein Logical, whether to generate representative protein FASTA
#'   files when `protein_sep` is supplied.
#'
#' @returns A list with class `ggsynteny_input_plan`.
#' @export
#'
#' @examples NULL
prepare_mcscanx_inputs <- function(design,
                                   outdir,
                                   protein_sep = NULL,
                                   prepare_protein = TRUE) {
  design <- validate_analysis_design(design)

  mcscanx_comparisons <- design$comparisons |>
    dplyr::filter(.data$analysis_tool == "mcscanx")

  if (nrow(mcscanx_comparisons) == 0) {
    stop("No comparisons in `design` use `analysis_tool = 'mcscanx'`.", call. = FALSE)
  }

  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }

  genomes_dir <- file.path(outdir, "genomes")
  comparisons_dir <- file.path(outdir, "comparisons")
  dir.create(genomes_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(comparisons_dir, recursive = TRUE, showWarnings = FALSE)

  genome_inputs <- design$genomes |>
    dplyr::mutate(
      genome_dir = file.path(genomes_dir, .data$genome_id),
      mcscanx_gff = file.path(.data$genome_dir, stringr::str_c(.data$genome_id, ".mcscanx.gff")),
      representative_protein = file.path(.data$genome_dir, stringr::str_c(.data$genome_id, ".representative.pep.fa")),
      protein_source = dplyr::if_else(
        !is.na(.data$protein_fasta) & !is.null(protein_sep) & prepare_protein,
        .data$representative_protein,
        .data$protein_fasta
      )
    )

  purrr::walk(genome_inputs$genome_dir, dir.create, recursive = TRUE, showWarnings = FALSE)

  purrr::pwalk(
    list(genome_inputs$annotation_path, genome_inputs$mcscanx_gff),
    function(annotation_path, mcscanx_gff) {
      if (!is.na(annotation_path) && nzchar(annotation_path)) {
        parse_gff3_to_mcscanx(gffpath = annotation_path, output = mcscanx_gff)
      }
    }
  )

  if (!is.null(protein_sep) && prepare_protein) {
    purrr::pwalk(
      list(genome_inputs$protein_fasta, genome_inputs$representative_protein),
      function(protein_fasta, representative_protein) {
        if (!is.na(protein_fasta) && nzchar(protein_fasta)) {
          get_representative_transcript(
            fasta = protein_fasta,
            sep = protein_sep,
            output = representative_protein
          )
        }
      }
    )
  }

  comparison_inputs <- mcscanx_comparisons |>
    dplyr::left_join(
      genome_inputs |>
        dplyr::select(
          query_genome = genome_id,
          query_mcscanx_gff = mcscanx_gff,
          query_protein = protein_source
        ),
      by = "query_genome"
    ) |>
    dplyr::left_join(
      genome_inputs |>
        dplyr::select(
          subject_genome = genome_id,
          subject_mcscanx_gff = mcscanx_gff,
          subject_protein = protein_source
        ),
      by = "subject_genome"
    ) |>
    dplyr::mutate(
      comparison_dir = file.path(comparisons_dir, .data$comparison_id),
      blast_output = file.path(.data$comparison_dir, stringr::str_c(.data$comparison_id, ".blast")),
      mcscanx_prefix = file.path(.data$comparison_dir, .data$comparison_id),
      combined_gff = file.path(.data$comparison_dir, stringr::str_c(.data$comparison_id, ".gff"))
    ) |>
    dplyr::arrange(.data$order)

  purrr::walk(comparison_inputs$comparison_dir, dir.create, recursive = TRUE, showWarnings = FALSE)

  purrr::pwalk(
    list(comparison_inputs$query_mcscanx_gff, comparison_inputs$subject_mcscanx_gff, comparison_inputs$combined_gff),
    function(query_gff, subject_gff, combined_gff) {
      if (!is.na(query_gff) && !is.na(subject_gff) && file.exists(query_gff) && file.exists(subject_gff)) {
        write_combined_text(c(query_gff, subject_gff), combined_gff)
      }
    }
  )

  structure(
    list(
      genomes = genome_inputs |>
        dplyr::select(
          genome_id,
          genome_label,
          order,
          annotation_path,
          protein_fasta,
          mcscanx_gff,
          representative_protein,
          protein_source
        ),
      comparisons = comparison_inputs |>
        dplyr::select(
          comparison_id,
          query_genome,
          subject_genome,
          order,
          query_protein,
          subject_protein,
          query_mcscanx_gff,
          subject_mcscanx_gff,
          combined_gff,
          blast_output,
          mcscanx_prefix,
          comparison_dir
        ),
      metadata = list(
        outdir = normalizePath(outdir, winslash = "/", mustWork = FALSE),
        analysis_tool = "mcscanx",
        protein_sep = protein_sep
      )
    ),
    class = "ggsynteny_input_plan"
  )
}

write_combined_text <- function(paths, output) {
  lines <- purrr::map(paths, readr::read_lines)
  readr::write_lines(unlist(lines, use.names = FALSE), output)
}
