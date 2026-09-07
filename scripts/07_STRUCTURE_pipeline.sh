#!/bin/bash
## ============================================================
## STRUCTURE pipeline: Bayesian clustering, K = 1-6, 10 replicates/K
##
## P. lividus - Fuencaliente CO2 vent gradient (La Palma)
##
## Run separately for the neutral SNP dataset and for the candidate
## SNPs under pH selection (as with PCA and FST), following the
## parallelized "one K per directory" approach: one directory per K,
## each running its 10 replicates in parallel with GNU parallel
## (adapted from command_PL_opt.sh, which only ran K = 6).
##
## Requires: the "structure" executable (v2.3.4), mainparams,
## extraparams, and GNU parallel installed.
##
## Usage:
##   1. Fill in the CONFIG section below (paths + NUMINDS/NUMLOCI
##      per dataset).
##   2. bash 05_STRUCTURE_pipeline.sh setup      # creates all K directories
##   3. Launch each dataset/K directory's run script (see note below)
##   4. bash 05_STRUCTURE_pipeline.sh collect     # gathers outputs into results/
## ============================================================

set -euo pipefail

## ------------------------------------------------------------
## CONFIG - edit before running
## ------------------------------------------------------------

STRUCTURE_BIN="./structure"          # path to the compiled structure executable
MAINPARAMS_TEMPLATE="./mainparams"   # template mainparams (neutral dataset settings)
EXTRAPARAMS="./extraparams"          # shared, unchanged across datasets/K

MIN_K=1
MAX_K=6          # see Materials and Methods: K = 1-6
NUMREPS=10        # replicates per K
THREADS=10        # parallel jobs (1 per replicate)

# One entry per dataset: label:infile:numinds:numloci
# NUMINDS/NUMLOCI must match each dataset's .structure file
DATASETS=(
  "neutral:populationsNeutro.structure:95:27310"
  "selection:populationsUS_RDA2_5.structure:95:232"
)

## ------------------------------------------------------------
## setup: create one directory per dataset x K, with all files
## structure needs, plus a run script for that K (parallelized
## over replicates with GNU parallel)
## ------------------------------------------------------------

setup() {
  for entry in "${DATASETS[@]}"; do
    IFS=':' read -r label infile numinds numloci <<< "$entry"

    for ((K=MIN_K; K<=MAX_K; K++)); do
      dir="${label}_K${K}"
      mkdir -p "$dir"

      cp "$STRUCTURE_BIN" "$dir/structure"
      cp "$EXTRAPARAMS" "$dir/extraparams"
      cp "$infile" "$dir/"

      # Build mainparams for this dataset (INFILE/NUMINDS/NUMLOCI substituted)
      sed \
        -e "s|^#define INFILE.*|#define INFILE  ${infile}  // (str) name of input data file|" \
        -e "s|^#define NUMINDS.*|#define NUMINDS    ${numinds}    // (int) number of diploid individuals in data file|" \
        -e "s|^#define NUMLOCI.*|#define NUMLOCI    ${numloci}    // (int) number of loci in data file|" \
        "$MAINPARAMS_TEMPLATE" > "$dir/mainparams"

      # Per-K run script: runs this K's 10 replicates in parallel
      cat > "$dir/run_K.sh" << INNER_EOF
#!/bin/bash
# Runs STRUCTURE for K=${K}, dataset "${label}", ${NUMREPS} replicates,
# ${THREADS} in parallel (1 thread per replicate). Each replicate gets
# its own random seed via -D \$RANDOM.
for ((REP=1; REP<=${NUMREPS}; REP++)); do
    echo "./structure -K ${K} -o output_k${K}.${label}\${REP} -D \$RANDOM"
done | parallel -j ${THREADS}
INNER_EOF
      chmod +x "$dir/run_K.sh"
    done
  done

  echo "Setup complete: one directory per dataset x K (e.g. neutral_K1 ... neutral_K6, selection_K1 ... selection_K6)."
  echo "Each directory has structure, mainparams, extraparams, the input file, and run_K.sh."
}

## ------------------------------------------------------------
## collect: after all runs finish, gather outputs per dataset
## into results/<dataset>/, ready to zip and upload to
## StructureSelector / CLUMPAK
## ------------------------------------------------------------

collect() {
  for entry in "${DATASETS[@]}"; do
    IFS=':' read -r label infile numinds numloci <<< "$entry"
    results_dir="results_${label}"
    mkdir -p "$results_dir"
    for ((K=MIN_K; K<=MAX_K; K++)); do
      mv "${label}_K${K}"/output_k${K}.${label}* "$results_dir/" 2>/dev/null || true
    done
    echo "Collected outputs for '${label}' into ${results_dir}/"
  done
}

## ------------------------------------------------------------
## main
## ------------------------------------------------------------

case "${1:-}" in
  setup)   setup ;;
  collect) collect ;;
  *) echo "Usage: $0 {setup|collect}"; exit 1 ;;
esac
