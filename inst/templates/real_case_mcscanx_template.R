library(ggSynteny)

config <- list(
  mcscanx_gff = "path/to/species_pair.gff",
  collinearity = "path/to/species_pair.collinearity",
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

data <- prepare_mcscanx_data(
  synteny = config$collinearity,
  gff = config$mcscanx_gff,
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

result <- validate_real_synteny_workflow(
  mcscanx_gff = config$mcscanx_gff,
  collinearity = config$collinearity,
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
