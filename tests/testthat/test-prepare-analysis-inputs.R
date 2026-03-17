test_that("prepare_mcscanx_inputs preserves genome and comparison order", {
  genomes <- register_genomes(
    genome_id = c("Ref", "Sp1", "Sp2"),
    annotation = c(
      extdata_path("toy_annotation.gff3"),
      extdata_path("toy_annotation.gff3"),
      extdata_path("toy_annotation.gff3")
    ),
    order = c(1, 2, 3)
  )

  comparisons <- define_comparisons(
    genomes = genomes,
    query = c("Ref", "Ref"),
    subject = c("Sp2", "Sp1"),
    order = c(2, 1),
    analysis_tool = "mcscanx"
  )

  design <- create_analysis_design(genomes, comparisons)
  outdir <- file.path(tempdir(), "ggsynteny-input-plan")

  plan <- prepare_mcscanx_inputs(
    design = design,
    outdir = outdir,
    prepare_protein = FALSE
  )

  expect_s3_class(plan, "ggsynteny_input_plan")
  expect_equal(plan$genomes$genome_id, c("Ref", "Sp1", "Sp2"))
  expect_equal(plan$comparisons$comparison_id, c("Ref_vs_Sp1", "Ref_vs_Sp2"))
  expect_equal(plan$comparisons$query_genome, c("Ref", "Ref"))
  expect_equal(plan$comparisons$subject_genome, c("Sp1", "Sp2"))
  expect_true(all(file.exists(plan$genomes$mcscanx_gff)))
  expect_true(all(file.exists(plan$comparisons$combined_gff)))
})

test_that("prepare_mcscanx_inputs creates combined comparison gff files", {
  genomes <- register_genomes(
    genome_id = c("A", "B"),
    annotation = c(extdata_path("toy_annotation.gff3"), extdata_path("toy_annotation.gff3")),
    order = c(1, 2)
  )

  comparisons <- define_comparisons(
    genomes = genomes,
    query = "A",
    subject = "B",
    analysis_tool = "mcscanx"
  )

  design <- create_analysis_design(genomes, comparisons)
  outdir <- file.path(tempdir(), "ggsynteny-input-plan-single")

  plan <- prepare_mcscanx_inputs(
    design = design,
    outdir = outdir,
    prepare_protein = FALSE
  )

  combined_lines <- readr::read_lines(plan$comparisons$combined_gff[[1]])
  expect_true(length(combined_lines) > 0)
  expect_true(any(stringr::str_detect(combined_lines, "A_gene1")))
  expect_true(any(stringr::str_detect(combined_lines, "B_gene1")))
})
