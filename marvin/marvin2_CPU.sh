#!/bin/bash
#SBATCH --job-name=lf
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
##SBATCH --time=24:00:00
#SBATCH --output=notebook-%j.out
#SBATCH --error=notebook-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=julia.vicens@upf.edu

echo "🔍 Estoy en el nodo: $(hostname)"
echo "🔁 Inicio del script: $(date)"

# Cargar módulo de Python necesario para Jupyter
module use /homes/aplic/noarch/modules/all
module load Python/3.11.5-GCCcore-13.2.0

# Añadir Julia a PATH
export PATH="/homes/users/jvicens/julia/julia-1.10.3/bin:$PATH"
echo "📦 Julia version:"
julia --version

# Comprobación de Jupyter
echo "📍 Ubicación de Jupyter:"
which jupyter
echo "📦 Jupyter version:"
jupyter --version

# Definir rutas del notebook
NOTEBOOK="../proves/Tod's/marvin/lessfriction/Tracking_c.ipynb"
OUTPUT="../proves/Tod's/marvin/lessfriction/Tracking_c_marvin.ipynb"

export JULIA_NUM_THREADS=64
# Ejecutar el notebook con kernel explícito
echo "🚀 Ejecutando notebook: $NOTEBOOK"
jupyter nbconvert --to notebook --execute "$NOTEBOOK" --output "$OUTPUT" \
  --ExecutePreprocessor.timeout=-1 \
  --ExecutePreprocessor.kernel_name=julia-1.10


echo "✅ Notebook ejecutado. Salida: $OUTPUT"
echo "🏁 Fin del script: $(date)"
