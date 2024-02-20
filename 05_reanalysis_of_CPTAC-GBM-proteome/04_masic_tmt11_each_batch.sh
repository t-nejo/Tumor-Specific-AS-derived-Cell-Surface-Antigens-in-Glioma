#!/bin/sh
#SBATCH --export=NONE
#SBATCH --mail-type=END,FAIL,REQEUE
#SBATCH --mail-user=takahide.nejo@ucsf.edu
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32gb
#SBATCH --time=24:00:00
#SBATCH --job-name=${JOB_NAME}
#SBATCH --output=${JOB_NAME}_%j.out
#SBATCH --error=${JOB_NAME}_%j.out
#SBATCH --gres=scratch:200G
#!/bin/sh

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
# input ; 
export D_IN_RAW=path_to_Proteome_Raw
export SAMPLE_LIST=path_to_sample_list.txt
export PARAM_FILE=path_to_param_file # masic_param_generate_ion_stats_tmt11.xml
# NOTE: this param file was generated from the default param file "MASICParameters.xml" with the following modifications. 
  # line 34: False -> True
  # line 35: value="0" -> value="16" 
# out ; 
export D_OUT=${D_CPTAC}/20210310_masic/out/${NAME1} ;

echo ' '
echo 'CHECK VARIABLES 2'
echo ' '
echo tmp_dir=${TMPDIR}
echo d_in_raw=${D_IN_RAW}
echo sample_list=${SAMPLE_LIST}
echo param_file=${PARAM_FILE}
echo d_out=${D_OUT}
echo ' '
echo 'PREP TMPDIR'
date
cd ${TMPDIR}
mkdir raw_files out
cp ${PARAM_FILE} ${TMPDIR}
echo ' '
echo 'START LOOP'
for i in \
$(seq 1 1 11) ; 
do
  echo 'SET VARIABLES'
  export NAME1=$(cat ${SAMPLE_LIST} | sed -n ${i}p | awk '{print $2}')
  echo ${i}
  echo name=${NAME1}
  echo  ' '
  echo 'CP FILES TO TMPDIR'
  echo ' '
  time cp \
  ${D_IN_RAW}/${NAME1}_f*.raw \
  ${TMPDIR}/raw_files
  echo ' '
  echo 'RUN MASIC'
  date ; 
  cd ${TMPDIR} ;

  time mono \
  path_to_masic/MASIC_Console.exe \
  /I:${TMPDIR}/raw_files/*.raw \
  /P:${TMPDIR}/${PARAM_FILE} \
  /O:${TMPDIR}/out \
  /L:${TMPDIR}/out/${NAME1}_log.txt \
  /SF:${TMPDIR}/out/${NAME1}_sf.xml ; 
  echo ' ' ;
done

echo 'CHECK'
cd ${TMPDIR}
ls -lh
echo ' ' 
tree 
echo 'MV OUTPUT'
date
if [ ! -d ${D_OUT} ]; then 
  mkdir -p ${D_OUT} ;
fi
cd ${TMPDIR}/out
mv * -t ${D_OUT}
echo ' '
echo 'MV OUTPUT DONE'
echo ' '
date
## End-of-job summary, if running as a job
[[ -n "$SLURM_JOB_ID" ]] && sstat --format="JobID,AveCPU,MaxRSS,MaxPages,MaxDiskRead,MaxDiskWrite" --allsteps --jobs="$SLURM_JOB_ID"
## end

