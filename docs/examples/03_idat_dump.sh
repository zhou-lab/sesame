## title: Peek at a raw IDAT -- idat-dump
set -eu
f=$SESAME_TEST_IDATS/EPICv2/206909630040_R03C01_Grn.idat.gz
## no `head` in the pipeline: it would SIGPIPE and trip `set -o pipefail`
sesame idat-dump --tsv "$f" > all_rows.tsv
awk 'NR<=5' all_rows.tsv > rows.tsv
## the page prints rows only: four numeric columns, no header
awk -F'\t' 'NR==1{exit !(NF==4 && $1 ~ /^[0-9]+$/)}' rows.tsv \
  || { echo "--tsv did not emit 4 numeric columns"; exit 1; }
## the summary form (no --tsv) is where the header lives
sesame idat-dump "$f" > summary.txt
grep -q 'addr' summary.txt || { echo "summary form lost its header"; exit 1; }
echo "idat-dump: $(wc -l < rows.tsv) rows, header only in the summary form"
