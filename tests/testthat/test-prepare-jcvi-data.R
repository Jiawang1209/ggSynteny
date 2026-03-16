test_that("prepare_jcvi_data returns a validated standard object", {
  dat <- prepare_jcvi_data(
    anchors = extdata_path("toy_jcvi.anchors"),
    bed_a = extdata_path("toy_jcvi_genomeA.bed"),
    bed_b = extdata_path("toy_jcvi_genomeB.bed"),
    genome1 = "GenomeA",
    genome2 = "GenomeB",
    gap_size = 50
  )

  expect_s3_class(dat, "ggsynteny_data")
  expect_equal(dat$metadata$source, "jcvi")
  expect_equal(nrow(dat$genomes), 2)
  expect_equal(nrow(dat$blocks), 2)
  expect_equal(sort(unique(dat$blocks$block_id)), c("0", "1"))
  expect_equal(sort(unique(dat$anchors$source)), "jcvi")
  expect_s3_class(validate_ggsynteny_data(dat), "ggsynteny_data")
})

test_that("plot_Synteny works on jcvi toy data", {
  dat <- prepare_jcvi_data(
    anchors = extdata_path("toy_jcvi.anchors"),
    bed_a = extdata_path("toy_jcvi_genomeA.bed"),
    bed_b = extdata_path("toy_jcvi_genomeB.bed"),
    genome1 = "GenomeA",
    genome2 = "GenomeB",
    gap_size = 50
  )

  p <- plot_Synteny(
    data = dat,
    min_pairs = 2,
    color_by = "block_id"
  )

  expect_s3_class(p, "ggplot")
})

test_that("micro-synteny plotting works on jcvi toy data without annotation features", {
  dat <- prepare_jcvi_data(
    anchors = extdata_path("toy_jcvi.anchors"),
    bed_a = extdata_path("toy_jcvi_genomeA.bed"),
    bed_b = extdata_path("toy_jcvi_genomeB.bed"),
    genome1 = "GenomeA",
    genome2 = "GenomeB",
    gap_size = 50
  )

  micro_dat <- prepare_micro_synteny_data(
    data = dat,
    block_id = "0",
    flank_genes = 1
  )

  expect_true("features" %in% names(micro_dat))
  expect_equal(nrow(micro_dat$features), 0)
  expect_true(all(c("feature_type", "plot_start", "plot_end") %in% names(micro_dat$features)))

  p <- plot_microSynteny(
    data = dat,
    block_id = "0",
    flank_genes = 1,
    show_features = TRUE
  )

  expect_s3_class(p, "ggplot")
})

test_that("jcvi workflow supports GFF3 annotation and transcript-aware micro-synteny", {
  res <- validate_jcvi_synteny_workflow(
    anchors = extdata_path("toy_jcvi.anchors"),
    bed_a = extdata_path("toy_jcvi_genomeA.bed"),
    bed_b = extdata_path("toy_jcvi_genomeB.bed"),
    genome1 = "GenomeA",
    genome2 = "GenomeB",
    gff3_genome1 = extdata_path("toy_annotation.gff3"),
    gff3_genome2 = extdata_path("toy_annotation.gff3"),
    transcript_mode = "canonical",
    min_pairs = 2,
    flank_genes = 1
  )

  expect_s3_class(res$data, "ggsynteny_data")
  expect_true("features" %in% names(res$data))
  expect_true("transcripts" %in% names(res$data))
  expect_s3_class(res$genome_plot, "ggplot")
  expect_s3_class(res$micro_plot, "ggplot")
})
