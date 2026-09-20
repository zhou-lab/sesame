## title: Read a .cg back -- attach-probe (the positional-storage card)
## A .cg holds no Probe_IDs; attach-probe joins the ordering back on. The page
## prints both the inferred-platform form and the named one.
set -eu
sesame attach-probe out/beta.cg > labelled.tsv
head -1 labelled.tsv | grep -q '^Probe_ID' || { echo "no Probe_ID header"; exit 1; }
grep -q '^cg00000029_TC21' labelled.tsv || { echo "cg00000029_TC21 missing"; exit 1; }

## naming the platform must give the same thing
sesame attach-probe --platform EPICv2 out/beta.cg > named.tsv
cmp -s labelled.tsv named.tsv || { echo "inferred and named disagree"; exit 1; }

## --all widens to every sample column. Note `head` in a pipeline would trip
## `set -o pipefail` with SIGPIPE, so the whole output goes to a file first.
sesame attach-probe --all out/beta.cg > all.tsv
awk -F'\t' 'NR==1{ if (NF!=3) exit 1; exit 0 }' all.tsv \
  || { echo "--all did not widen to 2 sample columns"; exit 1; }
echo "attach-probe: $(wc -l < labelled.tsv) rows, inferred == named"
