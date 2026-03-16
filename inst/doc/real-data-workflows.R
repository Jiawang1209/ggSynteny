library(ggSynteny)

mcscanx_gff <- system.file("extdata", "toy_mcscanx.gff", package = "ggSynteny")
collinearity <- system.file("extdata", "toy_mcscanx.collinearity", package = "ggSynteny")
annotation <- system.file("extdata", "toy_annotation.gff3", package = "ggSynteny")

res_mcscanx <- validate_real_synteny_workflow(
  mcscanx_gff = mcscanx_gff,
  collinearity = collinearity,
  genome1 = "GenomeA",
  genome2 = "GenomeB",
  gff3_genome1 = annotation,
  gff3_genome2 = annotation,
  transcript_mode = "canonical",
  min_pairs = 2,
  flank_genes = 1
)

names(res_mcscanx)
names(res_mcscanx$data)

jcvi_anchors <- system.file("extdata", "toy_jcvi.anchors", package = "ggSynteny")
bed_a <- system.file("extdata", "toy_jcvi_genomeA.bed", package = "ggSynteny")
bed_b <- system.file("extdata", "toy_jcvi_genomeB.bed", package = "ggSynteny")

res_jcvi <- validate_jcvi_synteny_workflow(
  anchors = jcvi_anchors,
  bed_a = bed_a,
  bed_b = bed_b,
  genome1 = "GenomeA",
  genome2 = "GenomeB",
  gff3_genome1 = annotation,
  gff3_genome2 = annotation,
  transcript_mode = "canonical",
  min_pairs = 2,
  flank_genes = 1
)
