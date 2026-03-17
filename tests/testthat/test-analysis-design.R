test_that("register_genomes preserves explicit order and metadata", {
  genomes <- register_genomes(
    genome_id = c("B", "A", "C"),
    annotation = c("B.gff3", "A.gff3", "C.gff3"),
    protein = c("B.pep.fa", "A.pep.fa", "C.pep.fa"),
    genome_label = c("Genome B", "Genome A", "Genome C"),
    order = c(2, 1, 3),
    group = c("grp1", "grp1", "grp2")
  )

  expect_equal(genomes$genome_id, c("A", "B", "C"))
  expect_equal(genomes$genome_label, c("Genome A", "Genome B", "Genome C"))
  expect_equal(genomes$order, c(1L, 2L, 3L))
})

test_that("define_comparisons preserves comparison direction and order", {
  genomes <- register_genomes(
    genome_id = c("A", "B", "C"),
    order = c(1, 2, 3)
  )

  comparisons <- define_comparisons(
    genomes = genomes,
    query = c("A", "A", "B"),
    subject = c("B", "C", "C"),
    order = c(2, 1, 3),
    analysis_tool = c("mcscanx", "jcvi", "mcscanx")
  )

  expect_equal(comparisons$query_genome, c("A", "A", "B"))
  expect_equal(comparisons$subject_genome, c("C", "B", "C"))
  expect_equal(comparisons$order, c(1L, 2L, 3L))
  expect_equal(comparisons$comparison_id, c("A_vs_C", "A_vs_B", "B_vs_C"))
})

test_that("create_analysis_design validates linked genome and comparison registries", {
  genomes <- register_genomes(
    genome_id = c("Ref", "Sp1", "Sp2"),
    order = c(1, 2, 3)
  )

  comparisons <- define_comparisons(
    genomes = genomes,
    query = c("Ref", "Ref"),
    subject = c("Sp1", "Sp2"),
    order = c(1, 2)
  )

  design <- create_analysis_design(genomes, comparisons)

  expect_s3_class(design, "ggsynteny_design")
  expect_equal(design$genomes$genome_id, c("Ref", "Sp1", "Sp2"))
  expect_equal(design$comparisons$query_genome, c("Ref", "Ref"))
  expect_equal(design$comparisons$subject_genome, c("Sp1", "Sp2"))
  expect_s3_class(validate_analysis_design(design), "ggsynteny_design")
})
