## title: Genotype SNPs -> VCF -- vcf
## The SNP table comes from the store: the page's --snp form is a PATH, with no
## store lookup behind it, which is why this runs without the flag.
set -eu
sesame vcf "$SESAME_TEST_IDATS/EPICv2/206909630040_R03C01" \
  --platform EPICv2 > geno.vcf
grep -q '^##fileformat=VCF' geno.vcf || { echo "not a VCF"; exit 1; }
n=$(grep -vc '^#' geno.vcf)
[ "$n" -gt 100 ] || { echo "only $n variant lines"; exit 1; }

sesame vcf "$SESAME_TEST_IDATS/EPICv2/206909630040_R03C01" \
  --platform EPICv2 --variants > informative.vcf
m=$(grep -vc '^#' informative.vcf)
[ "$m" -lt "$n" ] || { echo "--variants did not narrow the output"; exit 1; }
echo "vcf: $n sites, $m with --variants"
