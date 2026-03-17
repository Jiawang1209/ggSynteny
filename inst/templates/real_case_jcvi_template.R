library(ggSynteny)

config <- list(
  anchors = "path/to/species_pair.anchors",
  bed_a = "path/to/speciesA.bed",
  bed_b = "path/to/speciesB.bed",
  genome1 = "SpeciesA",
  genome2 = "SpeciesB",
  gff3_genome1 = "path/to/speciesA.gff3",
  gff3_genome2 = "path/to/speciesB.gff3",
  gene_id_map1 = NULL,
  gene_id_map2 = NULL,
  transcript_mode = "canonical",
  min_pairs = 5,
  flank_genes = 3
)

data <- prepare_jcvi_data(
  anchors = config$anchors,
  bed_a = config$bed_a,
  bed_b = config$bed_b,
  genome1 = config$genome1,
  genome2 = config$genome2
)

check_annotation_gene_id_match(
  data = data,
  gff3 = config$gff3_genome1,
  genome_id = config$genome1,
  gene_id_map = config$gene_id_map1
)

check_annotation_gene_id_match(
  data = data,
  gff3 = config$gff3_genome2,
  genome_id = config$genome2,
  gene_id_map = config$gene_id_map2
)

result <- validate_jcvi_synteny_workflow(
  anchors = config$anchors,
  bed_a = config$bed_a,
  bed_b = config$bed_b,
  genome1 = config$genome1,
  genome2 = config$genome2,
  gff3_genome1 = config$gff3_genome1,
  gff3_genome2 = config$gff3_genome2,
  gene_id_map1 = config$gene_id_map1,
  gene_id_map2 = config$gene_id_map2,
  transcript_mode = config$transcript_mode,
  min_pairs = config$min_pairs,
  flank_genes = config$flank_genes
)

result$genome_plot
result$micro_plot
