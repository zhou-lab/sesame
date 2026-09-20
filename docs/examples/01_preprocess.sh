## title: Preprocess -- IDAT to betas (the Preprocess card)
## The page's headline claim: QCDPB over a cohort in one process, writing one
## indexed .cg per output plus qc.tsv.
set -eu
idats=$SESAME_TEST_IDATS/EPICv2
sesame preprocess --out out/ \
  "$idats/206909630040_R03C01" "$idats/206909630042_R08C01"

for f in beta.cg intensity.cg pval.cg qc.tsv; do
  [ -s "out/$f" ] || { echo "missing out/$f"; exit 1; }
done
[ -s out/beta.cg.idx ] || { echo "beta.cg has no .idx of sample names"; exit 1; }

## the page says qc.tsv is 66 metrics beside `sample`
cols=$(head -1 out/qc.tsv | awk -F'\t' '{print NF}')
[ "$cols" = 67 ] || { echo "qc.tsv has $cols columns, expected 67"; exit 1; }
grep -q frac_dt out/qc.tsv || { echo "qc.tsv has no frac_dt column"; exit 1; }
echo "preprocess: 2 samples, $cols qc columns"
