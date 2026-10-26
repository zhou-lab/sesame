## title: Differential methylation -- dml (matrix form, and the .cg form)
## The model needs a real cohort; tests/run_dml.sh covers it against R. Here
## the documented matrix form, on six synthetic samples, which is what a reader
## can run without a cohort; then the page's .cg form with --platform EPICv2 on
## 01's two-sample out/beta.cg, which checks the Probe_IDs come from the store
## (two samples cannot fit the model, so the estimates are NA).
set -eu
awk 'BEGIN{OFS="\t";
  printf "Probe_ID"; for(s=1;s<=6;s++) printf "\tS%d", s; print "";
  srand(7);
  for(p=1;p<=200;p++){ printf "cg%07d", p;
    for(s=1;s<=6;s++){ b=(s<=3?0.3:0.7)+ (rand()-0.5)*0.05; printf "\t%.4f", b }
    print "" } }' > betas.tsv
printf 'sample\tgroup\n' > meta.tsv
for s in 1 2 3; do printf 'S%d\tA\n' "$s" >> meta.tsv; done
for s in 4 5 6; do printf 'S%d\tB\n' "$s" >> meta.tsv; done

sesame dml --betas betas.tsv --meta meta.tsv --formula '~ group' > dml.tsv
head -1 dml.tsv | grep -q 'Probe_ID' || { echo "no Probe_ID column"; exit 1; }
n=$(tail -n +2 dml.tsv | wc -l)
[ "$n" = 200 ] || { echo "dml wrote $n rows, expected 200"; exit 1; }
echo "dml: $n probes, $(head -1 dml.tsv | awk -F'\t' '{print NF}') columns"

printf 'sample\tgroup\n206909630040_R03C01\tA\n206909630042_R08C01\tB\n' \
  > meta.cg.tsv
sesame dml --betas out/beta.cg --platform EPICv2 \
  --meta meta.cg.tsv --formula '~ group' > dml.cg.tsv
n=$(tail -n +2 dml.cg.tsv | wc -l)
[ "$n" = 937690 ] || { echo "dml on the .cg wrote $n rows, expected 937690"; exit 1; }
sed -n '2p' dml.cg.tsv | cut -f1 | grep -q '^cg' ||
  { echo "dml on the .cg has no Probe_IDs"; exit 1; }
echo "dml .cg --platform EPICv2: $n probes"
