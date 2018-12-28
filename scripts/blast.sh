#/usr/bin/env bash

for j in al br cr es
do
    echo $j
    outfile=$j.blast.tab
    echo -e "qseqid\tsseqid\tqlen\tslen\tlength\tevalue\tppos\tqstart\tqend\tsstart\tsend\tsframe" > $outfile
    tblastn \
        -db blastdb/$j.fna \
        -outfmt "6 qseqid sseqid qlen slen length evalue ppos qstart qend sstart send sframe" \
        -num_threads 4 \
        -query yc.faa >> $outfile
done
