## title: Copy number -- cnv
set -eu
idats=$SESAME_TEST_IDATS/EPICv2
sesame preprocess --prep "" --raw-signal --output total_intensity \
  --out t/ "$idats/206909630040_R03C01"
sesame cnv --platform EPICv2 \
  t/total_intensity.cg segments.tsv bins.tsv
