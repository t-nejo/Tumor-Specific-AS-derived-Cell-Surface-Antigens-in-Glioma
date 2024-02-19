#!/bin/sh
#SBATCH --export=NONE
#SBATCH --mail-type=END,FAIL,REQEUE
#SBATCH --mail-user=takahide.nejo@ucsf.edu
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=200gb
#SBATCH --time=24:00:00
#SBATCH --job-name=${JOB_NAME}
#SBATCH --output=${JOB_NAME}_%j.out
#SBATCH --error=${JOB_NAME}_%j.out
#SBATCH --gres=scratch:300G

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
export MANIFEST=${manifest.tsv}
export D_REF=${path_to_flair_ref}
export GENOME=GRCh37.primary_assembly.genome_chr7.fa
export GTF=gencode.v29lift37.annotation_chr7_plus_PTPRZ1_alt_v7.gtf
export D_FQ=${path_to_guppy_call_fastq_output}
export FQ=${XXX}_all_merged.fastq.gz
export D_BED=${path_to_flair_correct_output}
export BED=${XXX}_all_flair_corrected_v7_all_corrected.bed
# output
export D_OUT1=${path_to_flair_collapse_output}
export ISOFORM=${XXX}_flair_collapsed.isoforms.fa
export OUT_PREFIX1=${XXX}_flair_collapsed
export D_OUT2=${path_to_flair_quantify_output}
export OUT_PREFIX2=${XXX}_flair_quantified
echo ' '
echo 'CHECK VARIABLES 2'
echo tmpdir=${TMPDIR}
echo manifest=${MANIFEST}
echo d_ref=${D_REF}
echo genome=${GENOME}
echo gtf=${GTF}
echo d_fq=${D_FQ}
echo fq=${FQ}
echo d_bed=${D_BED}
echo bed=${BED}
echo d_out1=${D_OUT1}
echo isoform=${ISOFORM}
echo out_prefix1=${OUT_PREFIX1}
echo d_out2=${D_OUT2}
echo out_prefix2=${OUT_PREFIX2}
echo ' '
echo 'cp files to tmpdir' ;
cd ${TMPDIR} ;
time cp ${D_REF}/${GENOME} ${D_REF}/${GTF} \
${D_REF}/${MANIFEST} \
${D_FQ}/*fastq.gz ${D_BED}/*_all_corrected.bed . 
echo ' '
echo 'CHECK'
tree -L 3 . ;
echo ' '
echo 'CONCATENATE FASTQ AND BED'
cd ${TMPDIR} ;
cat *fastq.gz > ${FQ}
cat *bed > ${BED}
echo 'CONCATENATE FASTQ AND BED DONE'
echo ' '
echo 'CHECK'
ls -lh ;
tree -L 3 . ;
echo ' '
echo 'module load'
date
module load CBI miniconda3
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
echo 'RUN FLAIR COLLAPSE'
cd ${TMPDIR}
mkdir output
date
time flair collapse \
-g ${GENOME} \
-q ${BED}  \
-r ${FQ} \
-t ${SLURM_NTASKS} \
-f ${GTF} \
-o output/${OUT_PREFIX1} \
--keep_intermediate \
--generate_map \
--check_splice \
--annotation_reliant generate \
--trust_ends
echo ' '
echo 'RUN FLAIR COLLAPSE DONE'
echo ' '
echo 'CHECK'
ls -lh output
tree -L 3 . 
echo ' '
echo 'MV OUTPUT'
if [ ! -d ${D_OUT1} ]; then 
  mkdir -p ${D_OUT1} ;
fi ;
cd ${TMPDIR}/output
cp ${ISOFORM} ${TMPDIR} ;
time mv ${OUT_PREFIX1}* -t ${D_OUT1}
echo ' '
echo 'MV OUTPUT DONE'
echo ' '
echo 'RUN FLAIR QUANTIFY'
cd ${TMPDIR}
mkdir output2
date
time flair quantify \
-r ${MANIFEST} \
-i ${ISOFORM} \
-o output2/${OUT_PREFIX2} \
-t ${SLURM_NTASKS} \
--generate_map \
--sample_id_only \
--trust_ends
echo ' '
echo 'RUN FLAIR QUANTIFY DONE'
echo ' '
echo 'CHECK'
ls -lh output
tree -L 3 . 
echo ' '
echo 'MV OUTPUT'
if [ ! -d ${D_OUT2} ]; then 
  mkdir -p ${D_OUT2} ;
fi ;
cd ${TMPDIR}/output2
time mv ${OUT_PREFIX2}* -t ${D_OUT2}
echo ' '
echo 'MV OUTPUT DONE'
echo ' '
echo 'CHECK'
cd ${TMPDIR}
ls -lh ; 
echo ' ' 
ls -lh output
echo ' '
ls -lh output2
echo ' '
tree -L 3 .
echo ' '
date
## End-of-job summary, if running as a job
[[ -n "$SLURM_JOB_ID" ]] && sstat --format="JobID,AveCPU,MaxRSS,MaxPages,MaxDiskRead,MaxDiskWrite" --allsteps --jobs="$SLURM_JOB_ID"
## end
