#!/usr/bin/env bash
echo hi
ls *_ITSx.summary.txt
(for txt in $(ls *_ITSx.summary.txt); do
    grep -P '\t' $txt  | tr -d ' ' |
        sed 's/.*file:/itsx/; s/.*ITSbyITSx:/its/;s/.*mainstrand:/its_plus/; ;s/.*tarystrand:/its_minus/; ;s/.*chimericbyITSx:/chimeric/; s/://' |
        csvtk transpose -tT
done;) | csvtk -tT concat

#{output.sum};
#        "seqkit stat -aT data/{x}.fa {x}_ITSx.{{full,LSU,SSU,5_8S}}.fasta > {output.sta}
