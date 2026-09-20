## title: Copy number -- cnv
## norun: needs a normal panel, which is not published (build with `make cnv-normals`)
set -eu
idats=$SESAME_TEST_IDATS/EPICv2
sesame preprocess --prep "" --raw-signal --output total_intensity \
  --out t/ "$idats/206909630040_R03C01"
sesame cnv --platform EPICv2 --normals EPICv2.cnvnormals.cg \
  t/total_intensity.cg segments.tsv bins.tsv
