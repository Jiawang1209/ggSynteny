# ggSynteny

`ggSynteny` is an R package for preparing MCScanX input and output files and
drawing genome-level synteny and local micro-synteny plots in R.

The package currently supports a standardized pairwise workflow:

1. Read an MCScanX `gff` file.
2. Read an MCScanX `.collinearity` file.
3. Convert them into standardized chromosome and block tables.
4. Optionally enrich the object with `GFF3` annotations.
5. Draw genome-level or transcript-aware micro-synteny plots.

The package also includes a first `jcvi` parser that converts `BED + anchors`
input into the same internal `ggsynteny_data` structure.

## Status

Implemented:

- Register multiple genomes and ordered pairwise comparisons
- Convert `GFF3` to MCScanX `gff` format
- Extract representative protein transcripts
- Run `BLASTP` and `MCScanX`
- Parse MCScanX and jcvi pairwise synteny inputs
- Build a standardized `ggsynteny_data` object
- Plot genome-level synteny views
- Plot transcript-aware micro-synteny with optional exon/CDS features
- Validate end-to-end workflows on bundled toy datasets

Planned next:

- Example datasets from real genomes
- More input adapters and multi-genome layouts
- Vignettes and broader real-data validation

## Installation

This repository is still under active development. For now, install locally:

```r
devtools::load_all(".")
```

## Minimal example

The package includes a synthetic MCScanX example in `inst/extdata`.

```r
library(ggSynteny)

gff_path <- system.file("extdata", "toy_mcscanx.gff", package = "ggSynteny")
collinearity_path <- system.file(
  "extdata",
  "toy_mcscanx.collinearity",
  package = "ggSynteny"
)

synteny_data <- prepare_mcscanx_data(
  synteny = collinearity_path,
  gff = gff_path,
  genome1 = "GenomeA",
  genome2 = "GenomeB",
  gap_size = 50
)

plot_Synteny(
  data = synteny_data,
  min_pairs = 2,
  color_by = "block_id",
  chr_color = "grey20",
  block_alpha = 0.45,
  genome_order = c("GenomeA", "GenomeB"),
  chr_order = list(
    GenomeA = c("A_chr2", "A_chr1"),
    GenomeB = c("B_chr2", "B_chr1")
  )
)
```

## jcvi example

The package includes a synthetic jcvi example in `inst/extdata`.

```r
jcvi_data <- prepare_jcvi_data(
  anchors = system.file("extdata", "toy_jcvi.anchors", package = "ggSynteny"),
  bed_a = system.file("extdata", "toy_jcvi_genomeA.bed", package = "ggSynteny"),
  bed_b = system.file("extdata", "toy_jcvi_genomeB.bed", package = "ggSynteny"),
  genome1 = "GenomeA",
  genome2 = "GenomeB",
  gap_size = 50
)

plot_Synteny(
  data = jcvi_data,
  min_pairs = 2,
  color_by = "block_id"
)
```

If you also have a full `GFF3` annotation, you can enrich the internal object
with strand and exon/CDS structure before plotting micro-synteny:

```r
annotation_path <- system.file("extdata", "toy_annotation.gff3", package = "ggSynteny")

synteny_data <- add_gff3_annotation(synteny_data, annotation_path, genome_id = "GenomeA")
synteny_data <- add_gff3_annotation(synteny_data, annotation_path, genome_id = "GenomeB")

plot_microSynteny(
  data = synteny_data,
  block_id = "0",
  flank_genes = 1,
  show_features = TRUE,
  show_gene_label = TRUE
)
```

## Real-data validation

When you have real MCScanX output and optional `GFF3` annotations, the fastest
way to validate the full workflow is:

```r
result <- validate_real_synteny_workflow(
  mcscanx_gff = "path/to/species_pair.gff",
  collinearity = "path/to/species_pair.collinearity",
  genome1 = "SpeciesA",
  genome2 = "SpeciesB",
  gff3_genome1 = "path/to/speciesA.gff3",
  gff3_genome2 = "path/to/speciesB.gff3",
  transcript_mode = "longest_cds",
  min_pairs = 5,
  flank_genes = 3
)

result$genome_plot
result$micro_plot
```

Required inputs:

- An MCScanX `gff` file matching the collinearity result
- An MCScanX `.collinearity` file
- Optional `GFF3` files if you want strand, transcript, exon, or CDS-aware micro-synteny

If the gene IDs in `GFF3` do not exactly match the MCScanX gene IDs, pass
`gene_id_map1` and/or `gene_id_map2` as named character vectors.

The same pattern now works for `jcvi` input:

```r
result <- validate_jcvi_synteny_workflow(
  anchors = "path/to/species_pair.anchors",
  bed_a = "path/to/speciesA.bed",
  bed_b = "path/to/speciesB.bed",
  genome1 = "SpeciesA",
  genome2 = "SpeciesB",
  gff3_genome1 = "path/to/speciesA.gff3",
  gff3_genome2 = "path/to/speciesB.gff3",
  transcript_mode = "canonical",
  min_pairs = 5,
  flank_genes = 3
)

result$genome_plot
result$micro_plot
```

Before attaching annotation on a real dataset, you can quickly inspect gene ID
compatibility with:

```r
check_annotation_gene_id_match(
  data = result$data,
  gff3 = "path/to/speciesA.gff3",
  genome_id = "SpeciesA"
)
```

Drop-in real-data script templates are bundled in:

- `inst/templates/real_case_mcscanx_template.R`
- `inst/templates/real_case_jcvi_template.R`

## Example data format

`toy_mcscanx.gff`

```text
A_chr1  A_gene1  1    100
A_chr1  A_gene2  160  260
B_chr1  B_gene1  10   110
B_chr1  B_gene2  180  300
```

`toy_mcscanx.collinearity`

```text
# synthetic example
0- 0:  A_gene1  B_gene1  1e-20
0- 1:  A_gene2  B_gene2  1e-25
```

## Main functions

- `prepare_mcscanx_data()`: parse MCScanX output into standardized tables
- `prepare_jcvi_data()`: parse jcvi `BED + anchors` output into standardized tables
- `register_genomes()`: register ordered genome inputs for multi-genome analysis
- `define_comparisons()`: define ordered pairwise comparisons
- `prepare_mcscanx_inputs()`: materialize MCScanX input files from an analysis design
- `plot_Synteny()`: draw a genome-level synteny plot
  Supports block filtering, genome order, chromosome order, and multiple color mappings.
- `plot_microSynteny()`: draw a block-level micro-synteny plot with local gene tracks
- `validate_real_synteny_workflow()`: run the full MCScanX validation path
- `validate_jcvi_synteny_workflow()`: run the full jcvi validation path
- `parse_gff3_to_mcscanx()`: convert `GFF3` to MCScanX input format
- `run_blastp_to_mcscanx()`: run BLASTP for MCScanX
- `run_MCScanX()`: run the MCScanX executable

## Development notes

- The package metadata and plotting API are still being refined.
- The current plotting API is designed around two genomes and standardized pairwise output.
- A vignette template for real-data runs is available in `vignettes/real-data-workflows.Rmd`.
