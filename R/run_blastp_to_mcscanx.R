run_blastp_to_mcscanx <- function(query, 
                                  db,
                                  output,
                                  blastdir,
                                  outfmt = 6,
                                  threads = 6,
                                  num_alignments = 5,
                                  evalue = 1e-10
                                  ){
  # blastdir = "/Users/liuyue/miniconda3/envs/blastp/bin"
  # db = "./Genome_data_2/Athaliana_167_protein.clean.fa"
  # query = "./Genome_data_2/Athaliana_167_protein.clean.fa"
  # output = "./Genome_data_2/Athaliana.blast"
  # outfmt = 6
  # threads = 6
  # num_alignments = 5
  # evalue = 1e-10
  
  # makeblastdb
  cmd1 <- sprintf("%s/makeblastdb -in %s -dbtype prot -title -parse_seqids -out %s", blastdir, db, db)
  
  system(cmd1)

  # blastp
  
  cmd2 <- sprintf(
    "%s/blastp -query %s -db %s -outfmt %d -out %s -num_threads %d -num_alignments %d -evalue %g",
    blastdir, query, db, outfmt, output, threads, num_alignments, evalue
  )
  
  system(cmd2)
  
}
