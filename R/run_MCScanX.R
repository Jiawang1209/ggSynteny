#' Run MCScanX synteny analysis
#'
#' @param prefix Character
#' The prefix of input file for MCScanX software
#' @param MCScanXdir Character
#' Path of MCScanX software
#'
#' @returns None
#' @export
#'
#' @examples NULL
#' 
run_MCScanX <- function(prefix,
                        MCScanXdir
                        ){
  cmd1 <- sprintf("%s/MCScanX %s", MCScanXdir, prefix)
  system(cmd1)
  
}
