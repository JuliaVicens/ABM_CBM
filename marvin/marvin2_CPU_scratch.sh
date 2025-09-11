#!/bin/bash
#SBATCH --job-name=lf
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
##SBATCH --time=24:00:00
#SBATCH --output=notebook-%j.out
#SBATCH --error=notebook-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=julia.vicens@upf.edu

echo "🔍 Nodo: $(hostname)"
echo "🕒 Inicio: $(date)"

# Módulos
module use /homes/aplic/noarch/modules/all
module load Python/3.11.5-GCCcore-13.2.0

# Julia (si no usas módulo Julia, PATH local)
export PATH="/homes/users/jvicens/julia/julia-1.10.3/bin:$PATH"
julia --version

# Jupyter check
which jupyter && jupyter --version

# ----------------------------
# 1) Definir SCRATCH y job dir
# ----------------------------
SCRATCH_BASE="/gpfs/scratch/lab_ojalvo/jvicens"
JOB_DIR="$SCRATCH_BASE/job_${SLURM_JOB_ID:-$$}"
mkdir -p "$JOB_DIR"

# ----------------------------
# 2) Notebook origen y nombres
# ----------------------------
SRC_NOTEBOOK="/homes/users/jvicens/ABM_CBM/proves/Tod's/attraction/prove/Tracking_c.ipynb"
RUN_NOTEBOOK="$JOB_DIR/Tracking_c.ipynb"                # copia que se ejecuta
OUT_NAME="motiles_few.ipynb"                   # nombre del ipynb de salida

# Copiar notebook al scratch
cp -f "$SRC_NOTEBOOK" "$RUN_NOTEBOOK"

# Entrar en el scratch (¡muy importante!)
cd "$JOB_DIR" || exit 1

# Hilos Julia coherentes con cpus-per-task
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-64}
echo "🧵 JULIA_NUM_THREADS=$JULIA_NUM_THREADS"

# ----------------------------
# 3) Ejecutar el notebook en scratch
#    Usa --output (nombre) + --output-dir (carpeta)
# ----------------------------
echo "🚀 Ejecutando notebook en: $RUN_NOTEBOOK"
jupyter nbconvert --to notebook --execute "$RUN_NOTEBOOK" \
  --output "$OUT_NAME" \
  --output-dir "$JOB_DIR" \
  --ExecutePreprocessor.timeout=-1 \
  --ExecutePreprocessor.kernel_name="julia-1.10"

NB_STATUS=$?

# ----------------------------
# 4) Copiar solo el resultado final a HOME
# ----------------------------

echo "✅ Copiado a: $FINAL_HOME_DIR/$OUT_NAME"
echo "🏁 Fin: $(date)"
exit $NB_STATUS
