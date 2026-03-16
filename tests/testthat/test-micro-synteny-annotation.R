test_that("GFF3 annotation adds strand, transcripts, and features", {
  dat <- prepare_mcscanx_data(
    synteny = extdata_path("toy_mcscanx.collinearity"),
    gff = extdata_path("toy_mcscanx.gff"),
    genome1 = "GenomeA",
    genome2 = "GenomeB",
    gap_size = 50
  )

  dat <- add_gff3_annotation(
    dat,
    gff3 = extdata_path("toy_annotation.gff3"),
    genome_id = "GenomeA",
    transcript_mode = "canonical"
  )
  dat <- add_gff3_annotation(
    dat,
    gff3 = extdata_path("toy_annotation.gff3"),
    genome_id = "GenomeB",
    transcript_mode = "canonical"
  )

  expect_true("features" %in% names(dat))
  expect_true("transcripts" %in% names(dat))
  expect_true(all(c("+", "-") %in% unique(stats::na.omit(dat$genes$strand))))
  expect_true(all(c("exon", "CDS") %in% unique(dat$features$feature_type)))

  selected <- dat$genes[, c("gene_id", "transcript_id")]
  expect_equal(selected$transcript_id[selected$gene_id == "A_gene1"], "A_tx1")
  expect_equal(selected$transcript_id[selected$gene_id == "B_gene1"], "B_tx1")
})

test_that("micro synteny preparation and plotting support transcript modes", {
  dat <- prepare_mcscanx_data(
    synteny = extdata_path("toy_mcscanx.collinearity"),
    gff = extdata_path("toy_mcscanx.gff"),
    genome1 = "GenomeA",
    genome2 = "GenomeB",
    gap_size = 50
  )

  dat <- add_gff3_annotation(dat, extdata_path("toy_annotation.gff3"), "GenomeA")
  dat <- add_gff3_annotation(dat, extdata_path("toy_annotation.gff3"), "GenomeB")

  micro_all <- prepare_micro_synteny_data(dat, block_id = "0", transcript_mode = "all")
  micro_custom <- prepare_micro_synteny_data(
    dat,
    block_id = "0",
    transcript_mode = "custom",
    transcript_ids = c(A_gene1 = "A_tx1", A_gene2 = "A_tx2", B_gene1 = "B_tx1", B_gene2 = "B_tx2")
  )

  expect_true(nrow(micro_all$features) >= nrow(micro_custom$features))
  expect_true(all(unique(micro_custom$features$transcript_id) %in% c("A_tx1", "A_tx2", "B_tx1", "B_tx2")))

  p <- plot_microSynteny(
    data = dat,
    block_id = "0",
    transcript_mode = "longest_cds",
    show_features = TRUE,
    show_gene_label = TRUE
  )

  expect_s3_class(p, "ggplot")
})
