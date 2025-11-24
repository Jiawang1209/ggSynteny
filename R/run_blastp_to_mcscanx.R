#' Run BLASTP for MCScanX synteny analysis
#'
#' @param query Character
#' Path to the protein FASTA file used as the BLASTP query.
#' @param db Character
#' Path to the protein FASTA file used to build the BLAST database.
#' @param output Character
#' Path to the output BLAST result file (typically `.blast`).
#' @param blastdir Character
#'  Directory containing BLAST+ executables
#' @param outfmt Character
#' BLAST output format (default: 6, tabular format).
#' @param threads Numeric
#' Number of threads for BLASTP (default: 6).
#' @param num_alignments Numeric
#' Maximum number of alignments to report per query (default: 5).
#' @param evalue Numeric
#' E-value threshold for BLASTP search (default: 1e-10).
#'
#' @returns a result of blastp
#' @export
#'
#' @examples NULL
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
