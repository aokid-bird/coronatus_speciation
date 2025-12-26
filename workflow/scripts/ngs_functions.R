# Minimal subset of ngs_functions needed for this pipeline

suppressPackageStartupMessages({
  library(tidyverse)
})

# Convert an ngsDist matrix into a NEXUS file
dist2nexus <- function(df,
                       labels = "no",
                       triangle = "both",
                       diagonal = TRUE,
                       dir,
                       dryrun = TRUE){
  ind <- df %>% dplyr::pull(ind) %>% length()
  taxlabels <- df %>% dplyr::pull(ind) %>% stringr::str_c(collapse = "\n")

  if(labels == "no"){
    df <- df %>% dplyr::select(!ind)
  }
  m <- df %>% apply(1, stringr::str_c, collapse = " ") %>% stringr::str_c(collapse = "\n")

  header <- "#nexus"
  taxablock <- stringr::str_c("BEGIN Taxa;",
                     glue::glue("DIMENSIONS ntax={ind};"),
                     "TAXLABELS",
                     taxlabels,
                     ";","END; [Taxa]",
                     sep = "\n")

  .diag <- ifelse(diagonal, "DIAGONAL", "NO DIAGONAL")

  distanceblock <- stringr::str_c("BEGIN Distances;",
                         glue::glue("DIMENSIONS ntax={ind};"),
                         glue::glue("FORMAT labels=LEFT { .diag } triangle={stringr::str_to_upper(triangle)};"),
                         "MATRIX",
                         m,
                         ";",
                         "END; [Distances]",
                         sep = "\n")

  nexus <- stringr::str_c(header, taxablock, distanceblock, sep = "\n\n")

  if(dryrun){
    cat(nexus)
  }else{
    readr::write_lines(nexus, dir)
  }
}

