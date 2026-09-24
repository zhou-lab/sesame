## title: Read a .cg back -- describe-probe (the positional-storage card)
## A .cg holds no Probe_IDs; describe-probe joins the ordering back on. The
## page prints both the inferred-platform form and the named one, then asks the
## same command where those probes are in the genome.
set -eu
sesame describe-probe out/beta.cg > labelled.tsv
head -1 labelled.tsv | grep -q '^Probe_ID' || { echo "no Probe_ID header"; exit 1; }
grep -q '^cg00000029_TC21' labelled.tsv || { echo "cg00000029_TC21 missing"; exit 1; }

## naming the platform must give the same thing
sesame describe-probe --platform EPICv2 out/beta.cg > named.tsv
cmp -s labelled.tsv named.tsv || { echo "inferred and named disagree"; exit 1; }

## --all widens to every sample column. Note `head` in a pipeline would trip
## `set -o pipefail` with SIGPIPE, so the whole output goes to a file first.
sesame describe-probe --all out/beta.cg > all.tsv
awk -F'\t' 'NR==1{ if (NF!=3) exit 1; exit 0 }' all.tsv \
  || { echo "--all did not widen to 2 sample columns"; exit 1; }
## --genome turns the same command around: probe IDs in, coordinates out, in
## the <chrm>_<beg1> spelling `yame rowsub -L` reads.
## awk, not `head`, because a pipeline that closes early trips pipefail
awk -F'\t' 'NR>1 && NR<=4 { print $1 }' labelled.tsv > probes.txt
sesame describe-probe --platform EPICv2 --genome hg38 probes.txt > coords.tsv
awk -F'\t' 'NF!=2 || $2 !~ /^chr[0-9XYM]+_[0-9]+$/ { exit 1 }' coords.tsv \
  || { echo "coordinate column is not <chrm>_<beg1>"; exit 1; }

echo "describe-probe: $(wc -l < labelled.tsv) rows, inferred == named, $(wc -l < coords.tsv) coordinates"
