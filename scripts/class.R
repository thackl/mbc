#!/usr/bin/env RScript
library(tidyverse)
library(gggenomes)

read_vsearch <- function(x, col_names=c('seq_id', 'lineage', 'strand', 'lineage_signif'), id="file_id"){
    read_tsv(x, col_names =  col_names, id=id)
}

tidy_lineage <- function(x, col){
  # str_match_all on NA returns matrix with one row of NA,NA
  # no match to "" gives empty matrix, which is what we need
  lng <- x |> pull({{col}}) |> replace_na("") 
  ss <- str_match_all(lng, "(\\b\\w)(?::|__)([^,;]+)")
  d <- tibble(.rows= length(ss))
  ii <- rep(seq_along(ss), map_int(ss, nrow))
  mm <- list_c(ss)
  kk <- factor(mm[,2])
  vv <- mm[,3]
  for (k in levels(kk)){
    d[[k]] <- NA
    d[[k]][ii[kk==k]] <- vv[kk==k]
  }

  bind_cols(x, d)
}

##----------------------------------------------------------------------------##
setwd("fungi-classify")
setwd("../../Frodo/n/")
setwd("ting-20230821/m")

vsearch_files <- list.files(".", "*sintax.tsv")
minimap_files <- list.files(".", "*.paf")

v0 <- read_vsearch(vsearch_files) |> select(-strand)
m0 <- read_feats(minimap_files) |> 
  group_by(file_id, seq_id) |> slice_head(n=1) |> ungroup() |> # highest scoring only
  select(file_id, seq_id, lineage=seq_id2)

d0 <- bind_rows(v0, m0) |>
  separate(file_id, into = c("sample_id", "region", "strategy"), sep="_" ) |>
  separate(seq_id, into = c("read_id", "region_info"), sep = "\\|", extra="merge")

boot_threshold <- .2
d1 <- d0 |>
  tidy_lineage(lineage) |> 
  mutate(across(7:13, function(.x){
    boot <- as.numeric(str_extract(.x, "[\\d.]+"))
    has_boot <- !is.na(boot)
    .x <- str_remove(.x, "\\(.*")
    .x <- ifelse(has_boot & boot < boot_threshold, "_LOW", .x)
    .x
    }))

# d k p c o f g s
#d1 |> filter(sample_id %in% c("frodo")) |> 
#d1 |> filter(str_detect(sample_id, "frodo(-|$)")) |> 
#d1 |> filter(sample_id %in% c("fa2303-bc01")) |> 
d1 |>
  ggplot() +
  # d k p c o f g s
  geom_histogram(aes(y=o, fill=paste(region, strategy)), stat="count", position=position_dodge(preserve="single")) +
  facet_wrap(~sample_id, nrow=1) +
  scale_x_sqrt() +
  scale_fill_brewer(palette="Paired") +
  theme_bw()


d1 |> 
  #filter(str_detect(sample_id, "frodo(-|$)")) |> 
  count(sample_id, region, strategy) |> 
  pivot_wider(names_from = c(region, strategy), values_from=n)

# ggsankey
# devtools::install_github("davidsjoberg/ggsankey")
library(ggalluvial)


da <- d1 |> filter(region=="ITS" & strategy == "vsearch-sintax.tsv") |> 
  count(k,p,c,o)

ggplot(data = da, aes(axis1 = k, axis2 = p, axis3 = c, y = n)) +
  scale_x_discrete(limits = c("k", "p", "c"), expand = c(.2, .05)) 
  geom_alluvium()

library(ggsankey)
df <- mtcars |> 
  make_long(cyl, vs, am, gear, carb)
df <- d1 |> filter(region=="ITS" & strategy == "vsearch-sintax.tsv") |> make_long(k,p,c) 
|> 
  arrange(node) |>
  mutate(next_node = factor(next_node, labels = unique(next_node)))


ggplot(df, aes(x = x, 
               next_x = next_x, 
               node = node, 
               next_node = next_node,
               fill = factor(node))) +
  geom_alluvial(flow.alpha = .6)
# ITSx & seq stats

ggplot() +
  geom_point(aes(x=read_id, y=o, shape=paste(region, strategy)), alpha=.3) +
  facet_wrap(~sample_id, scales="free_x") 

# read stats
itsx_files <- list.files(".", "*ITSx.summary.tsv")
i0 <- read_tsv(itsx_files, col_names = c("key", "n"), id="file_id")
i1 <- i0 |>
  mutate(file_id = str_remove(file_id, "_ITSx.summary.tsv")) |>
  pivot_wider(names_from="key", values_from="n")

seqs_files <- list.files(".", "*seqstats.tsv", full.names = T)
s0 <- read_tsv(seqs_files, id="file_id")
s1 <- s0 |> 
  mutate(
    region = str_extract(file, "[^.]+.fasta") |> str_remove(".fasta"),
    file_id = str_remove(basename(file_id), "_seqstats.tsv")) |> 
  select(file_id, region, n=num_seqs) |> 
  pivot_wider(names_from=region, values_from=n)

s2 <- left_join(s1, i1) |> relocate(file_id, itsx_in=itsx, itsx_any=its, its=full)

