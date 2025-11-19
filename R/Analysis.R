
####----step1:获取 代表性转录本----####
get_representative_transcript(fasta = "./Genome_data_2/Athaliana_167_protein.fa",
                              sep = ".",
                              output = "./Genome_data_2/Athaliana_167_protein.clean.fa")

get_representative_transcript(fasta = "./Genome_data_2/Osativa_323_v7.0.protein.fa",
                              sep = ".",
                              output = "./Genome_data_2/Osativa_323_v7.0.protein.clean.fa")

####----step2:解析gff3文件----####
parse_gff3_to_mcscanx(gffpath = "./Genome_data_2/Athaliana_167_gene.gff3",
                      output = "./Genome_data_2/Athaliana.gff")

parse_gff3_to_mcscanx(gffpath = "./Genome_data_2/Osativa_323_v7.0.gene.gff3",
                      output = "./Genome_data_2/Osativa.gff")

####----step3:获取染色体长度----####
parse_gff(gffpath = "./Genome_data_2/Athaliana_167_gene.gff3",
          output = "./Genome_data_2/Athaliana.Chr.info")

parse_gff(gffpath = "./Genome_data_2/Osativa_323_v7.0.gene.gff3",
          output = "./Genome_data_2/Osativa.Chr.info")

####----step4:进行blastp----####
# 要优先安装blast软件
# 自身 blastp
run_blastp_to_mcscanx(query = "./Genome_data_2/Athaliana_167_protein.clean.fa",
                      db = "./Genome_data_2/Athaliana_167_protein.clean.fa",
                      output = "./Genome_data_2/Athaliana.blast",
                      blastdir = "/Users/liuyue/miniconda3/envs/blastp/bin",
                      outfmt = 6,
                      threads = 6,
                      num_alignments = 5,
                      evalue = 1e-10)

# 多个物种的blastp
run_blastp_to_mcscanx(query = "./Genome_data_2/Athaliana_167_protein.clean.fa",
                      db = "./Genome_data_2/Osativa_323_v7.0.protein.clean.fa",
                      output = "./Genome_data_2/Athaliana_Osativa.blast",
                      blastdir = "/Users/liuyue/miniconda3/envs/blastp/bin",
                      outfmt = 6,
                      threads = 6,
                      num_alignments = 5,
                      evalue = 1e-10)

 # 在涉及到多个物种的时候，需要涉及的就是多个物种gff的合并
Athaliana_Osativa.blast <- dplyr::bind_rows(
  readr::read_delim(file = "./Genome_data_2/Athaliana.gff", col_names = F, delim = "\t") %>%
    dplyr::mutate(X1 = str_c("At", X1)),
  readr::read_delim(file = "./Genome_data_2/Osativa.gff", col_names = F, delim = "\t") %>%
    dplyr::mutate(X2 = str_remove(X2, pattern = ".MSUv7.0")) %>%
    dplyr::mutate(X1 = str_c("Os", X1))
)

write_delim(Athaliana_Osativa.blast,
            file = "./Genome_data_2/Athaliana_Osativa.gff",
            quote = "none",
            col_names = F,
            delim = "\t")


####----step5:MCScanX----####
# 对两个物种
run_MCScanX(prefix = "./Genome_data_2/Athaliana_Osativa",
            MCScanXdir = "/Users/liuyue/Desktop/Bioinformatics/MCScanX")


####----step6:解析共线性的结果----####
# 这里就需要解析一下共线性的结果了
synteny_block_out <- parse_Synteny_from_mcscanx(synteny = "Genome_data_2/Athaliana_Osativa.collinearity",
                                                gff = "Genome_data_2/Athaliana_Osativa.gff")

synteny_block_out[[1]]
synteny_block_out[[2]]

####----step7:可视化----####
# 首先优化一下染色体
# 染色体之间，我需要加一个gap
synteny_block_out[[2]] %>%
  dplyr::filter(stringr::str_detect(string = X1, pattern = "At")) %>%
  dplyr::slice(1:5) %>%
  dplyr::mutate()


# 我们单纯拿一个染色体出来试试
tmp_synteny <- synteny_block_out[[1]] %>% dplyr::slice(1:2)

tmp_synteny2 <- dplyr::bind_rows(
  tmp_synteny %>% dplyr::select(1:4) %>% purrr::set_names(c("block","Chr","Start","End")) %>% dplyr::mutate(species = "At"), 
  tmp_synteny %>% dplyr::select(1,5:7)%>% purrr::set_names(c("block","Chr","Start","End")) %>% dplyr::mutate(species = "Os")
)

tmp_chr <- synteny_block_out[[2]] %>% dplyr::slice(c(1,8)) %>%
  dplyr::mutate(Species = c("At", "Os")) %>%
  dplyr::mutate(Start = 1) %>%
  dplyr::rename(End = Length)

tmp_synteny2

tmp <- bind_rows(
  tmp_synteny2 %>% dplyr::filter(block == "0") %>% dplyr::slice(1) %>% dplyr::mutate(End = 3703672),
  tmp_synteny2 %>% dplyr::filter(block == "0") %>% dplyr::slice(1) %>% dplyr::mutate(Start = 3890811) ,
  tmp_synteny2 %>% dplyr::filter(block == "0") %>% dplyr::slice(2) %>% dplyr::mutate(End = 33000047),
  tmp_synteny2 %>% dplyr::filter(block == "0") %>% dplyr::slice(2) %>% dplyr::mutate(Start = 33460087),
  
  tmp_synteny2 %>% dplyr::filter(block == "1") %>% dplyr::slice(1) %>% dplyr::mutate(End = 22534717),
  tmp_synteny2 %>% dplyr::filter(block == "1") %>% dplyr::slice(1) %>% dplyr::mutate(Start = 22801073),
  tmp_synteny2 %>% dplyr::filter(block == "1") %>% dplyr::slice(2) %>% dplyr::mutate(End = 33000047),
  tmp_synteny2 %>% dplyr::filter(block == "1") %>% dplyr::slice(2) %>% dplyr::mutate(Start = 33361794)
  )

ggplot() + 
  geom_diagonal_wide(data = tmp, 
                     mapping = aes(x = Start, y = species, group = block, color = block, fill = block),
                     orientation = "y") + 
  ggnewscale::new_scale_fill() + 
  geom_gene_arrow(data = tmp_chr, mapping = aes(xmin = Start, xmax = End, y = Species, fill = Species),
                  arrowhead_width = grid::unit(0, "mm"),
                  arrowhead_height = grid::unit(0, "mm")) + 
  labs(x = "") + 
  theme_genes() + 
  theme(
    axis.text = element_text(color = "#000000", size = 15),
    axis.title = element_text(color = "#000000", size = 20),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  )



ggplot() + 
  geom_diagonal_wide(data = tmp, 
                     mapping = aes(x = Start, y = species, group = block, color = block, fill = block),
                     orientation = "y") + 
  ggnewscale::new_scale_fill() + 
  geom_gene_arrow(data = tmp_chr, mapping = aes(xmin = Start, xmax = End, y = Species, fill = Species),
                  arrowhead_width = grid::unit(0, "mm"),
                  arrowhead_height = grid::unit(0, "mm")) + 
  labs(x = "") + 
  theme_genes() + 
  theme(
    axis.text = element_text(color = "#000000", size = 15),
    axis.title = element_text(color = "#000000", size = 20),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  )


