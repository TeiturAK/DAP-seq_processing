#!/bin/bash

# Läs rad för rad från standard input (fil eller pipe)
while IFS= read -r r1_path || [ -n "$r1_path" ]; do
    # Hoppa över tomma rader
    [ -z "$r1_path" ] && continue

    # 1. Extrahera mappnamnet/provnamnet precis innan filnamnet
    # Exempel: /path/to/A_control/A_control_..._1.fq.gz -> A_control
    sample_id=$(basename "$(dirname "$r1_path")")

    # 2. Skapa söksträngen för Read 2 genom att ersätta _1.fq.gz med _2.fq.gz på slutet
    r2_path="${r1_path%_1.fq.gz}_2.fq.gz"

    # 3. Skriv ut CSV-raden
    echo "${sample_id},tf,${r1_path},${r2_path}"
done
