#!/bin/sh
#SBATCH --export=NONE
#SBATCH --mail-type=END,FAIL,REQEUE0
#SBATCH --mail-user=${YOUR_MAIL_ADRESS}
#SBATCH --nodes=1 
#SBATCH --ntasks=16
#SBATCH --mem=64gb
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
export FLOWCELL=FLO-MIN106
export KIT=SQK-LSK109
export BARCODE_KIT=EXP-NBD104
export NUM_CPU_THREADS=${SLURM_JOB_NUM_NODES}
export NUM_CPUS=${SLURM_NTASKS}
echo ' '
export i=${1} # 1, 2, ... 8
export ID=${ID} # sample.id
export D_FAST5=${path_to_fast5_directory} # directory where input fast5 files are stored.
export IN_FAST5=${FAM_XXX}_${ID}.fast5
# output
export D_OUT=${D_OUT} # directory to store the output 
echo ' '
echo 'CHECK VARIABLES 2'
echo flowcell=${FLOWCELL}
echo kit=${KIT}
echo barcode_kit=${BARCODE_KIT}
echo num_cpu_threads=${NUM_CPU_THREADS}
echo num_cpus=${NUM_CPUS}
echo i=${i}
echo id=${ID}
echo d_fast5=${D_FAST5}
echo in_fast5=${IN_FAST5}
echo d_out=${D_OUT}
echo ' '
echo 'CP FILES TO TMPDIR'
cd ${TMPDIR}
if [ ! -d ${ID} ] ; then 
  mkdir ${ID} ;
fi 
if [ ! -d input ] ; then 
  mkdir input ;
fi 
time cp -pR ${D_FAST5}/*/${IN_FAST5} ${TMPDIR}/input
echo ' '
echo 'RUN GUPPY_BASECALLER'
date
time guppy_basecaller \
--flowcell ${FLOWCELL} \
--kit ${KIT} \
--cpu_threads_per_caller ${NUM_CPU_THREADS} \
--num_callers ${NUM_CPUS} \
--compress_fastq \
--barcode_kits ${BARCODE_KIT} \
--trim_barcodes \
--qscore_filtering \
-i ${TMPDIR}/input \
-s ${TMPDIR}/${ID}
echo ' ' 
echo 'RUN GUPPY_BASECALLER DONE'
echo ' '
echo 'MV OUTPUT'
date
if [ ! -d ${D_OUT} ] ; then
  mkdir -p ${D_OUT} ; 
fi ; 
time mv ${TMPDIR}/${ID2} ${D_OUT}
echo ' '
echo 'MV OUTPUT DONE'
echo ' '
echo 'EDIT OUTPUT'
cd ${D_OUT}
for j in \
$(seq 1 1 8) ;
do 
  export BARCODE=barcode`printf %02d ${j}` ;
  if [ -e ${D_OUT}/${ID}/pass/${BARCODE}/*.fastq.gz ] ; then
    if [ ! -d ${D_OUT}/fastq_pass/${BARCODE} ] ; then 
      mkdir -p ${D_OUT}/fastq_pass/${BARCODE} ; 
    fi ; 
    mv ${D_OUT}/${ID}/pass/${BARCODE}/*.fastq.gz \
    ${D_OUT}/fastq_pass/${BARCODE}/FAM97239_pass_${BARCODE}_${ID}.fastq.gz ; 
  fi ;
done
echo ' '
if [ ! -e ${D_OUT}/sequencing_summary_merged.txt ] ; then
  cat ${D_OUT}/${ID}/sequencing_summary.txt > ${D_OUT}/sequencing_summary_merged.txt ;
else
  cat ${D_OUT}/${ID}/sequencing_summary.txt | tail -n +2 >> \
  ${D_OUT}/sequencing_summary_merged.txt ;
fi
if [ ! -d ${D_OUT}/done ] ; then 
  mkdir -p ${D_OUT}/done ; 
fi 
echo ' '
mv ${D_OUT}/${ID} -t ${D_OUT}/done ; 
echo ' '
echo 'EDIT OUTPUT DONE'
date
## End-of-job summary, if running as a job
[[ -n "$SLURM_JOB_ID" ]] && sstat --format="JobID,AveCPU,MaxRSS,MaxPages,MaxDiskRead,MaxDiskWrite" --allsteps --jobs="$SLURM_JOB_ID"
## end
