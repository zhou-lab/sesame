## title: Custom array -- an ordering by path, a mask by path (the Custom array tab)
## The tab's claim: a platform IS its ordering table, and --platform takes one
## by path as readily as by name. Built here from EPICv2's own ordering, so the
## result is checkable against the named-platform run.
set -eu
store=${YAME_DATA_HOME:-$HOME/.local/share/yame}
zcat "$store/EPICv2/EPICv2.ordering.tsv.gz" > MyArray.ordering.tsv

## the mask: one bit per ordering row, in the ordering's order, 1 = drop
awk -F'\t' 'NR>1{print ($1 ~ /^ctl_/) ? 1 : 0}' MyArray.ordering.tsv > m.txt
yame pack -f b m.txt MyArray.mask.cm
echo M_custom > names.txt
yame index -s names.txt MyArray.mask.cm

idats=$SESAME_TEST_IDATS/EPICv2
sesame preprocess --platform MyArray.ordering.tsv --mask MyArray.mask.cm \
  --prep QCDPB --out custom/ "$idats/206909630040_R03C01"
[ -s custom/beta.cg ] || { echo "custom-array preprocess wrote no betas"; exit 1; }

## and without a mask: prep C or "" work, but only with --output beta --
## the default --output carries pval and qc, which need a mask whatever the prep
sesame preprocess --platform MyArray.ordering.tsv --prep C --output beta \
  --out nomask/ "$idats/206909630040_R03C01"
[ -s nomask/beta.cg ] || { echo "mask-less run wrote no betas"; exit 1; }
if sesame preprocess --platform MyArray.ordering.tsv --prep C \
     --out shouldfail/ "$idats/206909630040_R03C01" 2>err.txt; then
  echo "expected the default --output to need a mask"; exit 1
fi
grep -q 'need a mask' err.txt || { echo "wrong error for a mask-less default run"; exit 1; }
echo "custom array: ordering by path, mask by path, and the documented refusal"
