test_that("prepare_mcscanx_data returns a validated standard object", {
  dat <- prepare_mcscanx_data(
    synteny = extdata_path("toy_mcscanx.collinearity"),
    gff = extdata_path("toy_mcscanx.gff"),
    genome1 = "GenomeA",
    genome2 = "GenomeB",
    gap_size = 50
  )

  expect_s3_class(dat, "ggsynteny_data")
  expect_true(all(c("genomes", "chromosomes", "genes", "anchors", "blocks", "metadata") %in% names(dat)))
  expect_equal(nrow(dat$genomes), 2)
  expect_equal(nrow(dat$blocks), 2)
  expect_equal(sort(unique(dat$blocks$block_id)), c("0", "1"))
  expect_equal(sort(unique(dat$anchors$genome_a)), "GenomeA")
  expect_equal(sort(unique(dat$anchors$genome_b)), "GenomeB")
  expect_s3_class(validate_ggsynteny_data(dat), "ggsynteny_data")
})

test_that("plot_Synteny returns a ggplot object on toy data", {
  dat <- prepare_mcscanx_data(
    synteny = extdata_path("toy_mcscanx.collinearity"),
    gff = extdata_path("toy_mcscanx.gff"),
    genome1 = "GenomeA",
    genome2 = "GenomeB",
    gap_size = 50
  )

  p <- plot_Synteny(
    data = dat,
    min_pairs = 2,
    color_by = "block_id",
    chr_order = list(
      GenomeA = c("A_chr2", "A_chr1"),
      GenomeB = c("B_chr2", "B_chr1")
    )
  )

  expect_s3_class(p, "ggplot")
})
