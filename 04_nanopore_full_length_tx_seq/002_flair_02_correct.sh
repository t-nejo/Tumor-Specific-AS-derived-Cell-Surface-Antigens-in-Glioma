#!/bin/sh
#SBATCH --export=NONE
#SBATCH --mail-type=END,FAIL,REQEUE
#SBATCH --mail-user=takahide.nejo@ucsf.edu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --job-name=${JOB_NAME}
#SBATCH --output=${JOB_NAME}_%j.out
#SBATCH --error=${JOB_NAME}_%j.out
#SBATCH --gres=scratch:200G

echo "Written by Takahide Nejo, MD, PhD"
echo "UC San Francisco, Neurological Surgery, Okada Lab"
echo ' '
date
echo ' '
echo 'SET TMPDIR'
if [[ -z "$TMPDIR" ]]; then
  if [[ -d /scratch ]]; then TMPDIR=/scratch/$USER; else TMPDIR=/tmp/$USER; fi
  mkdir -p "$TMPDIR"
  export TMPDIR
fi
echo ' '
echo 'CHECK VARIABLES 1'
echo ' '
echo slurm_job_id=${SLURM_JOB_ID}
echo slurm_job_num_nodes=${SLURM_JOB_NUM_NODES}
echo slurm_job_nodelist=${SLURM_JOB_NODELIST}
echo slurm_ntasks=${SLURM_NTASKS}
echo slurm_mem_per_nodes=${SLURM_MEM_PER_NODE}mb
echo num_cpu_threads=${NUM_CPU_THREADS}
echo num_cpus=${NUM_CPUS}
echo ' '
echo 'SET VARIABLES'
export i=${1}
export I=`printf %02d ${i}`
export D_REF=${path_to_flair_ref}
export GENOME=GRCh37.primary_assembly.genome_chr7.fa
export GTF=gencode.v29lift37.annotation_chr7_plus_PTPRZ1_alt_v7.gtf
export D_BED=${path_to_flair_align_output} # output of the previous step
export BED=${XXX}_flair_aligned.bed 
# output
export D_OUT=${path_to_flair_correct_output}
export OUT_PREFIX=${XXX}_flair_corrected # FAM97239_pass_barcode${I}_flair_corrected
echo ' '
echo 'CHECK VARIABLES 2'
echo tmpdir=${TMPDIR}
echo i=${i}
echo I=${I}
echo d_ref=${D_REF}
echo genome=${GENOME}
echo gtf=${GTF}
echo d_bed=${D_BED}
echo bed=${BED}
echo d_out=${D_OUT}
echo out_prefix=${OUT_PREFIX}
echo ' '
echo 'cp files to tmpdir'
cd ${TMPDIR}
time cp ${D_REF}/${GENOME} ${D_REF}/${GTF} ${D_BED}/${BED} .
echo ' '
echo 'CHECK'
tree -L 3 . 
echo ' '
echo 'module load'
date
module load CBI miniconda3 # depending on environments
echo ' '
echo 'CHECK VERSIONS'
echo 'conda --version'
conda --version
echo 'python --version'
python --version
echo ' '
echo 'ACTIVATE CONDA ENVIRONMENT'
conda activate flair
echo 'flair --version'
flair --version
echo ' '
echo 'RUN FLAIR CORRECT'
cd ${TMPDIR}
mkdir output
date
time flair correct \
-q ${BED} \
-g ${GENOME} \
-t ${SLURM_NTASKS} \
-f ${GTF} \
-o output/${OUT_PREFIX}
echo ' '
echo 'RUN FLAIR CORRECT DONE'
echo ' '
echo 'CHECK'
ls -lh output
tree -L 3 . 
echo ' '
echo 'MV OUTPUT'
if [ ! -d ${D_OUT} ]; then 
  mkdir -p ${D_OUT} ;
fi ;
cd ${TMPDIR}/output
time mv ${OUT_PREFIX}* ${D_OUT}
echo ' '
echo 'MV OUTPUT DONE'
echo ' '
date
## End-of-job summary, if running as a job
[[ -n "$SLURM_JOB_ID" ]] && sstat --format="JobID,AveCPU,MaxRSS,MaxPages,MaxDiskRead,MaxDiskWrite" --allsteps --jobs="$SLURM_JOB_ID"
## end

