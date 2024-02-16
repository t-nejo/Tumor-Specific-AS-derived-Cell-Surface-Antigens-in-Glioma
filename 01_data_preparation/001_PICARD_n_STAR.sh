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
export BAM_IN=${ID}_RNAseq.bam # filename of input BAM
export FQ1=${ID}_RNAseq_f1.fastq 
export FQ2=${ID}_RNAseq_f2.fastq
export D_REF=${REFERENCE_DIRECTORY} # directory where STAR-REF is located
export STAR_REF=star_index_Homo_sapiens_GRCh37_2.7.7a # STAR-REF directory
export BAM=${ID}_RNAseq_hg19_Aligned.sortedByCoord.out.bam # output BAM
export BAM_SORTED=${ID}_RNAseq_hg19_Aligned.sortedByCoord.out_sorted.bam # sorted BAM
echo ' '
echo 'CHECK VARIABLES 2'
echo ' '
echo id=${ID}
echo d_bam=${D_BAM}
echo bam_in=${BAM_IN}
echo fq1=${FQ1}
echo fq2=${FQ2}
echo d_ref=${D_REF}
echo star_ref=${STAR_REF}
echo bam=${BAM}
echo bam_sorted=${BAM_SORTED}
echo ' '
echo 'LOAD MODULES'
module load CBI star/2.7.7a samtools/1.10 picard/2.23.1
echo ' '
echo 'CHECK MODULES'
echo 'STAR --version'
STAR --version
echo 'samtools --version'
samtools --version
echo ' '
echo 'CP FILES TO TMPDIR'
cd ${TMPDIR}
mkdir input && cd input
time cp -pR ${D_REF}/${STAR_REF} .
time cp ${D_BAM}/${BAM_IN} .
echo ' '
echo 'CHECK'
cd ${TMPDIR}
ls -lh
echo ' ' 
tree -L 3 .
echo ' '
echo "RUN PICARD: convert bam to fastq"
date
cd ${TMPDIR}
time java -jar $PICARD SamToFastq \
I=input/${ID}_RNAseq.bam \
FASTQ=${FQ1} \
SECOND_END_FASTQ=${FQ2} \
VALIDATION_STRINGENCY=SILENT
echo ' '
echo "RUN PICARD DONE"
date
echo ' '
echo 'RUN STAR'
date 
cd ${TMPDIR}
time STAR \
--genomeDir input/${STAR_REF} \
--readFilesIn ${FQ1} ${FQ2} \
--runThreadN ${SLURM_NTASKS} \
--outFilterMultimapScoreRange 1 \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 10 \
--alignIntronMax 500000 \
--alignMatesGapMax 1000000 \
--sjdbScore 2 \
--alignSJDBoverhangMin 1 \
--genomeLoad NoSharedMemory \
--limitBAMsortRAM 80000000000 \
--outFilterMatchNminOverLread 0.33 \
--outFilterScoreMinOverLread 0.33 \
--sjdbOverhang 100 \
--outSAMstrandField intronMotif \
--outSAMattributes NH HI NM MD AS XS \
--limitSjdbInsertNsj 2000000 \
--outSAMunmapped None \
--outSAMtype BAM SortedByCoordinate \
--outSAMheaderHD @HD VN1.4 \
--twopassMode Basic \
--outSAMmultNmax 1 \
--outFileNamePrefix ${ID}_RNAseq_hg19_
echo ' '
echo 'RUN STAR DONE'
date
echo ' '
echo 'RUN SAMTOOLS SORT'
date
cd ${TMPDIR}
time samtools sort ${BAM} -o ${BAM_SORTED}
echo ' '
echo 'RUN SAMTOOLS SORT DONE'
date
echo ' '
echo 'RUN SAMTOOLS INDEX'
date
cd ${TMPDIR}
time samtools index -@ ${SLURM_NTASKS} \
${BAM_SORTED} \
${BAM_SORTED}.bai
echo ' '
echo 'RUN SAMTOOLS INDEX DONE'
date
echo ' '
echo 'DONE'
date
## End-of-job summary, if running as a job
[[ -n "$SLURM_JOB_ID" ]] && sstat --format="JobID,AveCPU,MaxRSS,MaxPages,MaxDiskRead,MaxDiskWrite" --allsteps --jobs="$SLURM_JOB_ID"
## end
