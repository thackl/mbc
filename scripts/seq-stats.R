#!/usr/bin/env -S Rscript --vanilla
# command line interface ------------------------------------------------------#
library(docopt)
library(tidyverse)
'seq-stats - aggregate seq-stats collected with seqkit stats

Usage:
  seq-stats.R <out.tsv> <seqkit-stats.tsv>...
' -> doc
opt <- docopt(doc)

#seq_files <- "/home/thackl/Research/Endophyte-complex/fungi-amplicon/fungi-classify/mbc/examples/n/seq-stats-detailed.tsv"
seq_files <- opt[["<seqkit-stats.tsv>"]]

s0 <- read_tsv(seq_files)
s1 <- s0 |> 
  mutate(
    sample = str_extract(file, "[^/]+$") |> 
      str_remove("(_ITSx.*|_filt)?.[^.]+(\\..z)?$"),
    step = str_extract(file, "[^.]+.fasta") |>
      str_remove(".fasta") |>
      str_replace("full", "its"),
    step = ifelse(str_detect(file, "_filt.fa"), "filtered", step) |> 
      replace_na("input")) |> 
  select(sample, step, n=num_seqs) |> 
  pivot_wider(names_from=step, values_from=n) |> 
  arrange(sample)
s1

write_tsv(s1, opt[["<out.tsv>"]])
