#!/bin/sh
#SBATCH --export=NONE
#SBATCH --mail-type=END,FAIL,REQEUE
#SBATCH --mail-user=takahide.nejo@ucsf.edu
#SBATCH --nodes=1
#SBATCH --ntasks=4
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
export D_MSGFP=path_to_msgfplus # program
export D_MZML=path_to_Proteome_PNNL_mzML
export FASTA=path_to_fasta_file # "combined" fasta
export CONFIG=path_to_config_file
export SAMPLE_LIST=path_to_sample_list.txt
# out ; 
export D_OUT=path_to_output_directory
echo ' '
echo 'CHECK VARIABLES 2'
echo ' '
echo tmp_dir=${TMPDIR}
echo d_msgfp=${D_MSGFP}
echo d_mzml=${D_MZML}
echo fasta=${FASTA}
echo config=${CONFIG}
echo sample_list=${SAMPLE_LIST}
echo d_out=${D_OUT}
echo ' '
echo 'PREP TMPDIR'
date
cd ${TMPDIR}
mkdir mzml input
cp ${FASTA} ${CONFIG} ${TMPDIR}/input
echo ' '
echo 'START LOOP'
for i in \
$(seq 1 1 11) ; 
do
  for j in \
  $(seq 1 1 24) ;
  do  
    echo 'SET VARIABLES'
    export j=${1} ;
    export i=${2} ; 
    export I=`printf %02d ${i}` ;
    export NAME1=$(cat ${SAMPLE_LIST} | sed -n ${j}p | awk '{print $1}')
    export NAME2=$(cat ${SAMPLE_LIST} | sed -n ${j}p | awk '{print $2}')
    export FN_OUT=$(echo ${NAME2}_f${I}.mzML | sed -e 's/.mzML//g') 
    echo ' '
    echo 'CHECK VARIABLES'
    echo i=${i}
    echo I=${I}
    echo j=${j}
    echo name1=${NAME1}
    echo name2=${NAME2}
    echo fn_out=${FN_OUT}
    echo ' '
    echo 'CP FILES TO TMPDIR'
    echo ' ' 
    time gunzip -c ${D_MZML}/${NAME1}/${NAME2}_f${I}.mzML.gz > \
    ${TMPDIR}/mzml/${NAME2}_f${I}.mzML
    echo ' '
    echo 'RUN MSGF+'
    date 
    cd ${TMPDIR}
    echo ' '
    time java -Xmx3500M -jar ${D_MSGFP}/MSGFPlus.jar \
    -s mzml/${NAME2}_f${I}.mzML \
    -d input/${FASTA} \
    -conf input/${CONFIG} \
    -o ${FN_OUT}.mzid
    echo ' '
    date
    echo 'RUN MSGF+ DONE'
    echo ' '
    echo 'RUN MzIDToTsv' 
    date
    time java -cp ${D_MSGFP}/MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv \
    -i ${FN_OUT}.mzid \
    -showQValue 1 \
    -showDecoy 1 \
    -showFormula 1 \
    -unroll 0 
    echo ' '
    date
    echo 'RUN MzIDToTsv DONE'
    echo ' ' ;
  done 
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
fi ;
mv *.tsv *.mzid -t ${D_OUT}
echo ' '
echo 'MV OUTPUT DONE'
echo ' '
date
## End-of-job summary, if running as a job
[[ -n "$SLURM_JOB_ID" ]] && sstat --format="JobID,AveCPU,MaxRSS,MaxPages,MaxDiskRead,MaxDiskWrite" --allsteps --jobs="$SLURM_JOB_ID"
## end
