#!/bin/bash
# One task per sector (array index -> SECTORS[i]); each task multithreads over
# every mode combo in MODES_FILE. Submit: sbatch slurm/run_sector.slurm
# Subset: sbatch --array=4,5 slurm/run_sector.slurm

#SBATCH --job-name=qscan
#SBATCH --array=0-5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=REPLACE_ME     # Julia threads; Della's AMD EPYC nodes have 192 cores/node
#SBATCH --mem-per-cpu=REPLACE_ME       # e.g. 4G
#SBATCH --time=REPLACE_ME              # e.g. 08:00:00
#SBATCH --output=slurm/logs/%x_%A_%a.out
#SBATCH --error=slurm/logs/%x_%A_%a.err
#SBATCH --mail-type=REPLACE_ME         # e.g. fail, or delete this line
#SBATCH --mail-user=REPLACE_ME         # <netid>@princeton.edu, or delete this line

set -euo pipefail
mkdir -p slurm/logs

module purge
module load julia/1.12.1
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# We parallelise with Julia threads over mode combos; nested OpenMP/BLAS threads
# would only oversubscribe. Also silences libgomp's warning about an inherited
# non-numeric OMP_NUM_THREADS.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

MODES_FILE="$PROJECT_DIR/slurm/modes.txt"
OUTDIR="/scratch/gpfs/REPLACE_ME/qscan_${SLURM_ARRAY_JOB_ID}"   # <ResearchGroup>/<netid>; not backed up

SAME_FREQ=true            # true: omega1 = omega2; false: independent (omega1,omega2) grid
W1_MIN=REPLACE_ME         # e.g. 0.05
QNM_FRAC=REPLACE_ME       # w_max = QNM_FRAC * Re(omega_QNM(l_final, n=0)), e.g. 0.8
N1=REPLACE_ME              # points in w1

W2_MIN=REPLACE_ME         # only used if SAME_FREQ=false
N2=REPLACE_ME               # only used if SAME_FREQ=false

LAMBDA1=4.0
LAMBDA2=5.0

SECTORS=(eee eeo eoe eoo ooe ooo)
SECTOR="${SECTORS[$SLURM_ARRAY_TASK_ID]}"

cd "$PROJECT_DIR"
if [ "$SAME_FREQ" = true ]; then
    julia --project=. --threads="$SLURM_CPUS_PER_TASK" jl/scan_modes.jl "$SECTOR" "$MODES_FILE" "$OUTDIR" true \
        "$W1_MIN" "$QNM_FRAC" "$N1" "$LAMBDA1" "$LAMBDA2"
else
    julia --project=. --threads="$SLURM_CPUS_PER_TASK" jl/scan_modes.jl "$SECTOR" "$MODES_FILE" "$OUTDIR" false \
        "$W1_MIN" "$QNM_FRAC" "$N1" "$W2_MIN" "$N2" "$LAMBDA1" "$LAMBDA2"
fi
