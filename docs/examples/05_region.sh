## title: Region & gene views -- region (the cinderplot input)
## The page pipes this into cinderplot; the part sesame owns is the long-form
## TSV, which is what this checks.
set -eu
sesame region --gene DNMT3A --pad 2000 \
  --betas out/beta.cg --platform EPICv2 > dnmt3a.tsv
[ "$(wc -l < dnmt3a.tsv)" -gt 1 ] || { echo "no probes in DNMT3A"; exit 1; }

sesame region chr7:27,090,000-27,210,000 \
  --betas out/beta.cg --platform EPICv2 > hoxa.tsv
[ "$(wc -l < hoxa.tsv)" -gt 1 ] || { echo "no probes in the HOXA window"; exit 1; }
echo "region: DNMT3A $(($(wc -l < dnmt3a.tsv)-1)) rows, HOXA $(($(wc -l < hoxa.tsv)-1)) rows"
