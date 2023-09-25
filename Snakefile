# relative to workflow dir
wf_dir=config["wf_dir"]
configfile: f"{wf_dir}/config.yaml"

import pandas as pd
import glob
import requests
import re
import os
from os.path import exists

scripts=os.path.join(wf_dir, "scripts")

dbs=config["db_dir"]
if not os.path.isabs(dbs):
    dbs=os.path.join(wf_dir, dbs)

# read samples
def read_samples():
    samples_tsv = config["samples"]
    print(samples_tsv)
    df = pd.read_csv(samples_tsv, sep="\t")
    files = df["#file"].values.tolist()
    # adjust paths relative to sample tsv
    paths = []
    for file in files:
        if not os.path.isabs(file):
           file=os.path.join(os.path.dirname(samples_tsv), file)
        paths.append(file)

    # get ids, remove extension with optional compression
    files = [os.path.basename(file) for file in files]
    ids = [re.sub("\.[^.]+(\..z)?$", "", file) for file in files]
    df["id"] = ids
    df["path"] = paths
    return(df.to_dict("list"))

SAMPLES = read_samples()
print(SAMPLES)
import sys

X = SAMPLES["id"]
print("Sample table:\n",  SAMPLES)
IDMAP = dict(zip(SAMPLES["id"], SAMPLES["path"]))
MIMAP = dict(zip(SAMPLES["id"], SAMPLES["min_length"]))
MAMAP = dict(zip(SAMPLES["id"], SAMPLES["max_length"]))
print(MIMAP)

## Pseudo-rules
rule all:
    input: expand("vsearch/{x}_{region}_vsearch-sintax.tsv", x=X, region=["ITS" ,"SSU"]),
           expand("minimap2/{x}_{region}_minimap2.paf", x=X, region=["ITS" ,"SSU"]),
           "seq-stats.tsv"
           

rule stats:
    input:
        SAMPLES["path"],
        expand("prefiltered/{x}_filt.fa", x=X),
        expand("itsx/{x}_ITSx.{region}.fasta", x=X, region=["full", "SSU", "LSU", "5_8S"])
    output:
        det="seq-stats-detailed.tsv",
        sum="seq-stats.tsv"
    shell:
        "seqkit stat -aT {input} > {output.det};"
        "{scripts}/seq-stats.R {output.sum} {output.det}"
           
## analysis -----------------------------------------------------------------##

rule prefilter:
    input: lambda wildcards: IDMAP[wildcards.x]
    output:
        fa="prefiltered/{x}_filt.fa"
    params:
        m = lambda wildcards: MIMAP[wildcards.x],
        M = lambda wildcards: MAMAP[wildcards.x]
    shell:
        "seqkit fq2fa {input} | "
        "seqkit seq -m {params.m} -M {params.M} -o {output}"

rule extract_marker_itsx:
    input: "prefiltered/{x}_filt.fa"
    output: multiext("itsx/{x}_ITSx", ".full.fasta", ".SSU.fasta", ".LSU.fasta", ".5_8S.fasta", ".summary.txt")
            # ".5_8S.fasta", ".chimeric.fasta", ".ITS1.fasta", ".ITS2.fasta", 
    params:
        pre="itsx/{x}_ITSx"
    threads: workflow.cores
    shell:
        "ITSx {config[itsx]} -i {input} -o {params.pre} --cpu {threads};"

# the ITSx summary is off for large files, might be due to chunking, ...
# rule summarize_marker_itsx:
#     input:
#         sum="{x}_ITSx.summary.txt",
#         fas=expand("{{x}}_ITSx.{region}.fasta", region=["full", "SSU", "LSU", "5_8S"])
#     output:
#         sum="{x}_ITSx.summary.tsv",
#         sta="{x}_seqstats.tsv"
#     shell:
#         "(grep -P '\t' {input.sum} | tr -d ' ' |"
#         "sed 's/.*file:/itsx/; s/.*ITSbyITSx:/its/;s/.*mainstrand:/its_plus/; ;s/.*tarystrand:/its_minus/; ;s/.*chimericbyITSx:/chimeric/; s/://'"
#         ")> {output.sum};"
#         "seqkit stat -aT {input.fas} > {output.sta}"
        

rule classify_its_vsearch:
    input: "itsx/{x}_ITSx.full.fasta"
    output: "vsearch/{x}_ITS_vsearch-sintax.tsv"
    shell:
        "vsearch --sintax {input} --db {dbs}/its-unite-sintax.udb --tabbedout {output}"

rule classify_ssu_vsearch:
    input: "itsx/{x}_ITSx.SSU.fasta"
    output: "vsearch/{x}_SSU_vsearch-sintax.tsv"
    shell:
        "vsearch --sintax {input} --db {dbs}/ssu-silva-sintax.udb --tabbedout {output}"

rule classify_its_minimap2:
    input: "itsx/{x}_ITSx.full.fasta"
    output: "minimap2/{x}_ITS_minimap2.paf"
    shell:
        "minimap2 -cx map-ont -n 5 --secondary=no {dbs}/its-unite-sintax.mmi {input} > {output}"

rule classify_ssu_minimap2:
    input: "itsx/{x}_ITSx.SSU.fasta"
    output: "minimap2/{x}_SSU_minimap2.paf"
    shell:
        "minimap2 -cx map-ont -n 5 --secondary=no {dbs}/ssu-silva-sintax.fa {input} > {output}"

    
