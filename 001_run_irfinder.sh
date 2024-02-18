#!/bin/sh
#SBATCH --export=NONE
#SBATCH --mail-type=END,FAIL,REQEUE
#SBATCH --mail-user=${YOUR_MAIL_ADRESS}
#SBATCH --nodes=1 
#SBATCH --ntasks=8
#SBATCH --mem=128gb
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
echo ' '
export ID=${ID} # sample id
export D_BAM=${BAM_FILE_DIRECTORY} # directory where input BAM files are located
export BAM=${ID}_RNAseq_hg19_Aligned.sortedByCoord.out.bam
export BAI=${BAM}.bai
export D_IRF_REF=path/to/irfinder
export IRF_REF=Human-hg19-release75
export D_IRF_OUT=path/to/irfinder_output
echo ' '
echo 'CHECK VARIABLES 2'
echo ' '
echo id=${ID}
echo d_bam=${D_BAM}
echo bam=${BAM}
echo bai=${BAI}
echo d_irf_ref=${D_IRF_REF}
echo irf_ref=${IRF_REF}
echo d_irf_out=${D_IRF_OUT}
echo tmpdir=${TMPDIR}
echo ' '
echo "CP FILES TO TMPDIR"
date
time cp ${D_BAM}/${BAM} ${D_BAM}/${BAI} ${TMPDIR}
echo ' '
export IRF_REF_IN_TMPDIR=${TMPDIR}/${IRF_REF} ;
if [ ! -d ${IRF_REF_IN_TMPDIR} ]; then
  time cp -pR ${D_IRF_REF}/${IRF_REF} ${IRF_REF_IN_TMPDIR} ;
fi
echo ' '
cd ${TMPDIR}
if [ -e irf ]; then
  rm -rf irf
fi
mkdir irf
echo ' '
echo 'RUN IRFINDER'
date
cd ${TMPDIR}
time IRFinder -m BAM \
-r ${IRF_REF_IN_TMPDIR} \
-d irf \
${BAM}
echo 'RUN IRFINDER DONE'
date
echo ' '
echo 'CHECK'
ls -lh
echo ' '
echo 'CP FILES'
cd ${TMPDIR} 
cp irf/IRFinder-IR-nondir.txt ${D_IRF_OUT}/${ID}_IRFinder-IR-nondir.txt
ceho ' '
echo 'DONE'
date
## End-of-job summary, if running as a job
[[ -n "$SLURM_JOB_ID" ]] && sstat --format="JobID,AveCPU,MaxRSS,MaxPages,MaxDiskRead,MaxDiskWrite" --allsteps --jobs="$SLURM_JOB_ID"
## end
