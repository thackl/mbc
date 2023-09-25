#!/usr/bin/env -S Rscript --vanilla
# command line interface ------------------------------------------------------#
library(docopt)
library(tidyverse)
'seq-stats - aggregate seq-stats collected with seqkit stats

Usage:
  seq-stats.R <out.tsv> <seqkit-stats.tsv>...
' -> doc
opt <- docopt(doc)
#read_tsv("/home/thackl/Research/Endophyte-complex/fungi-amplicon/fungi-classify/Frodo/n/seq-stats-detailed.tsv")
seq_files <- opt[["<seqkit-stats.tsv>"]]
s0 <- read_tsv(seq_files)
s1 <- s0 |> 
  mutate(
    region = str_extract(file, "[^.]+.fasta") |>
      str_remove(".fasta") |>
      str_replace("full", "its") |> 
      str_replace("^", "itsx_") |>  
      replace_na("total"),
    sample = str_remove(basename(file), "\\.fa$|_ITSx.*")) |> 
  select(sample, region, n=num_seqs) |> 
  pivot_wider(names_from=region, values_from=n) |> 
  arrange(sample)

write_tsv(s1, opt[["<out.tsv>"]])