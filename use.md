# ggSynteny Use Guide

这份文档按“先跑通，再细看 API”的方式整理当前包的全部用法。

## 目录

1. 快速开始
2. 包的核心对象
3. 多基因组设计层
4. MCScanX 工作流
5. jcvi 工作流
6. 注释增强和 gene ID 检查
7. Genome-level 绘图
8. Micro-synteny 绘图
9. 一步跑完整流程
10. 输入预处理
11. 外部工具封装
12. 坐标和数据层辅助函数
13. 模板脚本
14. 导出函数索引
15. 推荐使用顺序

## 1. 快速开始

### 1.1 加载包

开发阶段建议直接：

```r
devtools::load_all(".")
library(ggSynteny)
```

### 1.2 toy 数据路径

```r
mcscanx_gff <- system.file("extdata", "toy_mcscanx.gff", package = "ggSynteny")
collinearity <- system.file("extdata", "toy_mcscanx.collinearity", package = "ggSynteny")
toy_gff3 <- system.file("extdata", "toy_annotation.gff3", package = "ggSynteny")

jcvi_anchors <- system.file("extdata", "toy_jcvi.anchors", package = "ggSynteny")
jcvi_bed_a <- system.file("extdata", "toy_jcvi_genomeA.bed", package = "ggSynteny")
jcvi_bed_b <- system.file("extdata", "toy_jcvi_genomeB.bed", package = "ggSynteny")
```

### 1.3 最快跑通 MCScanX

```r
dat <- prepare_mcscanx_data(
  synteny = collinearity,
  gff = mcscanx_gff,
  genome1 = "GenomeA",
  genome2 = "GenomeB",
  gap_size = 50
)

plot_Synteny(dat, min_pairs = 2)
```

### 1.4 最快跑通带注释的 micro-synteny

```r
dat <- add_gff3_annotation(dat, toy_gff3, genome_id = "GenomeA")
dat <- add_gff3_annotation(dat, toy_gff3, genome_id = "GenomeB")

plot_microSynteny(
  data = dat,
  block_id = "0",
  flank_genes = 1,
  show_features = TRUE
)
```

## 2. 包的核心对象

`ggSynteny` 的正式主入口都会返回一个统一对象：`ggsynteny_data`。

典型结构：

```r
names(dat)
```

通常包含：

- `genomes`
- `chromosomes`
- `genes`
- `anchors`
- `blocks`
- `metadata`

如果接入了 `GFF3` 注释，还可能包含：

- `features`
- `transcripts`

验证对象是否合规：

```r
validate_ggsynteny_data(dat)
```

## 3. 多基因组设计层

### 3.1 注册多个 genome

```r
genomes <- register_genomes(
  genome_id = c("Ref", "Sp1", "Sp2"),
  annotation = c("Ref.gff3", "Sp1.gff3", "Sp2.gff3"),
  protein = c("Ref.pep.fa", "Sp1.pep.fa", "Sp2.pep.fa"),
  order = c(1, 2, 3)
)
```

### 3.2 定义比较顺序

```r
comparisons <- define_comparisons(
  genomes = genomes,
  query = c("Ref", "Ref"),
  subject = c("Sp1", "Sp2"),
  order = c(1, 2),
  analysis_tool = "mcscanx"
)
```

### 3.3 组合成分析设计

```r
design <- create_analysis_design(genomes, comparisons)
validate_analysis_design(design)
```

### 3.4 根据 design 批量准备 MCScanX 输入

```r
plan <- prepare_mcscanx_inputs(
  design = design,
  outdir = "analysis_inputs",
  prepare_protein = FALSE
)

plan$genomes
plan$comparisons
```

## 4. MCScanX 工作流

### 3.1 正式主入口

```r
dat_mcscanx <- prepare_mcscanx_data(
  synteny = collinearity,
  gff = mcscanx_gff,
  genome1 = "GenomeA",
  genome2 = "GenomeB",
  gap_size = 50
)
```

### 3.2 底层读取函数

```r
genes_tbl <- read_mcscanx_gff(mcscanx_gff)
legacy <- parse_Synteny_from_mcscanx(collinearity, mcscanx_gff)
```

说明：

- `prepare_mcscanx_data()`：正式主入口，返回 `ggsynteny_data`
- `read_mcscanx_gff()`：只解析 MCScanX `gff`
- `parse_Synteny_from_mcscanx()`：兼容旧接口，返回 legacy list

## 5. jcvi 工作流

### 4.1 正式主入口

```r
dat_jcvi <- prepare_jcvi_data(
  anchors = jcvi_anchors,
  bed_a = jcvi_bed_a,
  bed_b = jcvi_bed_b,
  genome1 = "GenomeA",
  genome2 = "GenomeB",
  gap_size = 50
)
```

### 4.2 底层读取函数

```r
bed_tbl <- read_jcvi_bed(jcvi_bed_a, genome = "GenomeA")

legacy <- parse_Synteny_from_jcvi(
  anchors = jcvi_anchors,
  bed_a = jcvi_bed_a,
  bed_b = jcvi_bed_b,
  genome1 = "GenomeA",
  genome2 = "GenomeB"
)
```

说明：

- `prepare_jcvi_data()`：正式主入口，返回 `ggsynteny_data`
- `read_jcvi_bed()`：只读取 jcvi `BED`
- `parse_Synteny_from_jcvi()`：兼容旧接口

## 6. 注释增强和 gene ID 检查

### 5.1 读取 GFF3 注释

```r
anno <- read_gff3_annotation(toy_gff3)
names(anno)
```

返回：

- `genes`
- `transcripts`
- `features`

### 5.2 检查 gene ID 是否匹配

接真实数据前，建议先做这一步：

```r
match_info <- check_annotation_gene_id_match(
  data = dat_mcscanx,
  gff3 = toy_gff3,
  genome_id = "GenomeA"
)

match_info$match_rate
match_info$unmatched_data_gene_ids
match_info$unmatched_annotation_gene_ids
```

### 5.3 把 GFF3 接入标准对象

```r
dat_annot <- add_gff3_annotation(
  data = dat_mcscanx,
  gff3 = toy_gff3,
  genome_id = "GenomeA",
  transcript_mode = "canonical"
)

dat_annot <- add_gff3_annotation(
  data = dat_annot,
  gff3 = toy_gff3,
  genome_id = "GenomeB",
  transcript_mode = "canonical"
)
```

### 5.4 gene ID 不一致时的映射

```r
gene_id_map <- c(
  "gff_gene_1" = "mcscanx_gene_1",
  "gff_gene_2" = "mcscanx_gene_2"
)

dat_annot <- add_gff3_annotation(
  data = dat_mcscanx,
  gff3 = "speciesA.gff3",
  genome_id = "SpeciesA",
  gene_id_map = gene_id_map
)
```

### 5.5 transcript_mode 说明

当前支持：

- `"longest_transcript"`
- `"longest_cds"`
- `"canonical"`
- `"all"`
- `"custom"`

## 7. Genome-level 绘图

### 6.1 最小示例

```r
p1 <- plot_Synteny(
  data = dat_mcscanx,
  min_pairs = 2
)
```

### 6.2 常用参数示例

```r
p2 <- plot_Synteny(
  data = dat_mcscanx,
  min_pairs = 2,
  min_block_width = 50,
  color_by = "block_id",
  chr_color = "grey20",
  block_alpha = 0.45,
  genome_order = c("GenomeA", "GenomeB"),
  chr_order = list(
    GenomeA = c("A_chr2", "A_chr1"),
    GenomeB = c("B_chr2", "B_chr1")
  ),
  show_chr_label = TRUE,
  show_genome_label = TRUE
)
```

### 6.3 `plot_Synteny()` 常用参数

- `min_pairs`：最少 anchor 数过滤
- `min_block_width`：最小 block 宽度过滤
- `color_by`：颜色映射方式
- `genome_order`：上下 genome 顺序
- `chr_order`：每个 genome 的染色体顺序
- `chr_color`：染色体条带颜色
- `block_alpha`：block 连线透明度
- `show_chr_label`：是否显示染色体标签
- `show_genome_label`：是否显示 genome 标签

`color_by` 目前支持：

- `"n_anchors"`
- `"block_id"`
- `"genome_pair"`

## 8. Micro-synteny 绘图

### 7.1 先准备局部 block 数据

```r
micro_dat <- prepare_micro_synteny_data(
  data = dat_annot,
  block_id = "0",
  flank_genes = 1,
  transcript_mode = "canonical"
)
```

返回：

- `genomes`
- `genes`
- `anchors`
- `features`
- `block`
- `block_regions`

### 7.2 最小示例

```r
pm1 <- plot_microSynteny(
  data = dat_annot,
  block_id = "0"
)
```

### 7.3 常用参数示例

```r
pm2 <- plot_microSynteny(
  data = dat_annot,
  block_id = "0",
  flank_genes = 1,
  color_by = "block_id",
  gene_color_by = "anchor_status",
  show_features = TRUE,
  feature_types = c("exon", "CDS"),
  transcript_mode = "canonical",
  show_gene_label = TRUE,
  highlight_block = TRUE,
  block_fill = "#fee8c8"
)
```

### 7.4 自定义 transcript

```r
pm3 <- plot_microSynteny(
  data = dat_annot,
  block_id = "0",
  transcript_mode = "custom",
  transcript_ids = c(
    A_gene1 = "A_tx1",
    A_gene2 = "A_tx2",
    B_gene1 = "B_tx1",
    B_gene2 = "B_tx2"
  )
)
```

### 7.5 `plot_microSynteny()` 常用参数

- `block_id`：要画的 block
- `flank_genes`：两侧扩展的基因数
- `color_by`：anchor 连线颜色映射
- `gene_color_by`：gene 填充颜色映射
- `show_features`：是否显示 exon/CDS
- `feature_types`：显示哪些 feature 类型
- `transcript_mode`：转录本选择策略
- `show_gene_label`：是否标 anchor gene 标签
- `highlight_block`：是否高亮 block 区域
- `block_fill`：block 高亮颜色

`gene_color_by` 目前支持：

- `"anchor_status"`
- `"genome"`

## 9. 一步跑完整流程

### 8.1 MCScanX 完整验证

```r
res_mcscanx <- validate_real_synteny_workflow(
  mcscanx_gff = mcscanx_gff,
  collinearity = collinearity,
  genome1 = "GenomeA",
  genome2 = "GenomeB",
  gff3_genome1 = toy_gff3,
  gff3_genome2 = toy_gff3,
  transcript_mode = "canonical",
  min_pairs = 2,
  flank_genes = 1
)

res_mcscanx$data
res_mcscanx$genome_plot
res_mcscanx$micro_plot
```

### 8.2 jcvi 完整验证

```r
res_jcvi <- validate_jcvi_synteny_workflow(
  anchors = jcvi_anchors,
  bed_a = jcvi_bed_a,
  bed_b = jcvi_bed_b,
  genome1 = "GenomeA",
  genome2 = "GenomeB",
  gff3_genome1 = toy_gff3,
  gff3_genome2 = toy_gff3,
  transcript_mode = "canonical",
  min_pairs = 2,
  flank_genes = 1
)

res_jcvi$data
res_jcvi$genome_plot
res_jcvi$micro_plot
```

适用场景：

- 想最快验证一套输入是否能完整跑通
- 想直接拿到 `data + genome_plot + micro_plot`
- 想对真实数据先做一次端到端检查

## 10. 输入预处理

### 9.1 GFF3 转 MCScanX gff

```r
parse_gff3_to_mcscanx(
  gffpath = "species.gff3",
  output = "species.mcscanx.gff"
)
```

### 9.2 从 GFF3 提取染色体长度

```r
parse_gff(
  gffpath = "species.gff3",
  output = "chr_length.txt"
)
```

### 9.3 提取代表转录本

按最长蛋白序列保留一个 transcript：

```r
get_representative_transcript(
  fasta = "species.pep.fa",
  sep = ".",
  output = "species.pep.longest.fa"
)
```

## 11. 外部工具封装

### 10.1 运行 BLASTP

```r
run_blastp_to_mcscanx(
  query = "speciesA.pep.longest.fa",
  db = "speciesB.pep.longest.fa",
  output = "speciesA_vs_speciesB.blast",
  blastdir = "/path/to/blast/bin",
  outfmt = 6,
  threads = 8,
  num_alignments = 5,
  evalue = 1e-10
)
```

### 10.2 运行 MCScanX

```r
run_MCScanX(
  prefix = "species_pair",
  MCScanXdir = "/path/to/MCScanX"
)
```

这里的 `prefix` 指的是 MCScanX 运行所需那组同前缀文件。

## 12. 坐标和数据层辅助函数

这些函数一般不需要在常规使用时单独调用，但在你想自己组织数据时会有用。

```r
chr_layout <- compute_chr_layout(dat_mcscanx$chromosomes, gap_size = 50)
genes_linear <- map_genes_to_linear(dat_mcscanx$genes, chr_layout)
blocks_linear <- compute_block_layout(dat_mcscanx$blocks, chr_layout)
blocks_from_anchors <- build_blocks_from_anchors(dat_mcscanx$anchors)
```

说明：

- `compute_chr_layout()`：生成线性染色体布局
- `map_genes_to_linear()`：gene 到线性坐标
- `compute_block_layout()`：block 到线性坐标
- `build_blocks_from_anchors()`：由 anchors 汇总 block

## 13. 模板脚本

真实数据模板已经放在：

- `inst/templates/real_case_mcscanx_template.R`
- `inst/templates/real_case_jcvi_template.R`

通常只需要改：

- 文件路径
- `genome1` / `genome2`
- `gene_id_map1` / `gene_id_map2`
- `transcript_mode`
- `min_pairs`
- `flank_genes`

## 14. 导出函数索引

### 13.1 主工作流

- `register_genomes()`
- `define_comparisons()`
- `create_analysis_design()`
- `validate_analysis_design()`
- `prepare_mcscanx_inputs()`
- `prepare_mcscanx_data()`
- `prepare_jcvi_data()`
- `validate_real_synteny_workflow()`
- `validate_jcvi_synteny_workflow()`

### 13.2 读取和解析

- `read_mcscanx_gff()`
- `read_jcvi_bed()`
- `read_gff3_annotation()`
- `parse_Synteny_from_mcscanx()`
- `parse_Synteny_from_jcvi()`

### 13.3 注释和检查

- `add_gff3_annotation()`
- `check_annotation_gene_id_match()`
- `validate_ggsynteny_data()`

### 13.4 绘图

- `plot_Synteny()`
- `plot_microSynteny()`
- `prepare_micro_synteny_data()`

### 13.5 数据层辅助

- `build_blocks_from_anchors()`
- `compute_chr_layout()`
- `compute_block_layout()`
- `map_genes_to_linear()`

### 13.6 输入预处理

- `parse_gff3_to_mcscanx()`
- `parse_gff()`
- `get_representative_transcript()`

### 13.7 外部工具封装

- `run_blastp_to_mcscanx()`
- `run_MCScanX()`

### 13.8 其他导出

- `%>%`

## 15. 推荐使用顺序

### 14.1 常规分析

1. 准备 `MCScanX` 或 `jcvi` 输入
2. `register_genomes()`
3. `define_comparisons()`
4. `create_analysis_design()`
5. `prepare_mcscanx_inputs()` 或后续其他输入准备函数
6. `prepare_*_data()` 或读取分析结果
7. `check_annotation_gene_id_match()`
8. `add_gff3_annotation()`
9. `plot_Synteny()`
10. `plot_microSynteny()`

### 14.2 最快跑通

如果你只想最快拿到结果，直接用：

- `validate_real_synteny_workflow()`
- `validate_jcvi_synteny_workflow()`

### 14.3 接真实数据时最容易出问题的地方

- gene ID 对不上
- 染色体名对不上
- `GFF3` 不是同一版本基因组
- `GFF3` 缺 transcript/exon/CDS
- block anchor 太少，导致 micro-synteny 区域过 sparse
