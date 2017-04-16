#!/bin/bash -l
#$ -N QCathena
#$ -j y
#$ -m e
#$ -cwd
#$ -l athena=true
#$ -M ashley.doane@gmail.com
#$ -pe smp 4
#$ -l h_vmem=6G
#$ -R y
#$ -o /home/asd2007/joblogs

########### SETTINGS ###########
ATHENA=1
##############################

#cd /athena/elementolab/scratch/asd2007/Projects/DataSets/atacData/atacGiorgio

path=$1 #path to all the Samples
gtf_path=$2 #Number indicating the reference gtf [1 for Human 2 for Mouse 3 for other]

# Uses job array for each Sample in the folder
file=$(ls ${path} | tail -n +${SGE_TASK_ID}| head -1)

#change directory to node directory
#cd $TMPDIR
#cd ${file}
echo "File : $file Read from $path\n"

#Obtain the name of the Sample
Sample=$(basename "$file")
Sample=${Sample%%.*}



PICARDROOT="/home/asd2007/Tools/picard/build/libs"
export $PICARDROOT

cd ${path}
cd ${Sample}

ATHENA=1
gtf_path=1

if [ $ATHENA = 1 ]
then
    REFDIR="/athena/elementolab/scratch/asd2007/Reference"
    ANNOTDIR="/athena/elementolab/scratch/asd2007/Reference"
    #cp /home/asd2007/Scripts/picardmetrics.athena.conf ./picardmetrics.conf
else
    REFDIR="/zenodotus/dat01/melnick_bcell_scratch/asd2007/Reference"
    ANNOTDIR="/zenodotus/dat01/melnick_bcell_scratch/asd2007/Reference"
fi

if [ $gtf_path = 1 ]
then
    #gtf="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/mm10_UCSC_ref.gtf"
	REF="${REFDIR}/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex"
	REFbt2="${REFDIR}/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index"
  # REFbt="/home/asd20i07/dat02/asd2007/Reference/Homo_sapiens/UCSC/mm10/Sequence/BowtieIndex/genome"
  BLACK="${REFDIR}/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
  chrsz="/athena/elementolab/scratch/asd2007/Reference/hg19.chrom.sizes"
  cp $chrsz $PWD/hg19.chrom.sizes
  RG="hg19"
  SPEC="hs"
  REFGen="/athena/elementolab/scratch/asd2007/bin/bcbio/genomes/Hsapiens/hg19/seq/"
fi

 #   rm picardmetrics.conf

#    picardmetrics run -f "/home/asd2007/picardmetrics.athena.conf"  -o ${Sample}/QCmetrics/filtered ${Sample}/${Sample}.sorted.nodup.noM.bam
#cp ${Sample}/QCmetrics/filtered/*.EstimateLibraryComplexity.log ${Sample}/QCmetrics/${Sample}.picardcomplexity.qc
#cp ${Sample}/QCmetrics/filtered/*.EstimateLibraryComplexity.log ${Sample}/${Sample}.picardcomplexity.qc


export PICARD="/home/asd2007/Tools/picard/build/libs/picard.jar"

export PATH="/home/asd2007/Tools/picard/build/libs:$PATH"

JAVA_HOME=/home/akv3001/jdk1.8.0_05
#alias picard="java -jar $PICARD"

alias picard="java -Xms500m -Xmx6G -jar $PICARD"


ALIGNED_BAM="${PWD}/${Sample}/${Sample}.bam"
WORKDIR=$PWD
OUTDIR="qc"
OUTPREFIX=$Sample
INPREFIX=$Sample
GENOME='hg19' # This is the only genome that currently works

SAMPLE="$Sample"
# Annotation files

#ANNOTDIR="/home/asd2007/melnick_bcell_scratch/asd2007/Reference"
#ANNOTDIR="/athena/elementolab/scratch/asd2007/Reference"

DNASE_BED="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz"
BLACKLIST_BED="${ANNOTDIR}/${GENOME}/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
TSS_BED="${ANNOTDIR}/${GENOME}/hg19_RefSeq_stranded.bed.gz"
REF_FASTA="${ANNOTDIR}/${GENOME}/encodeHg19Male.fa"
PROM="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_prom_p2.bed.gz"
ENH="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_enh_p2.bed.gz"
REG2MAP="${ANNOTDIR}/${GENOME}/dnase_avgs_reg2map_p10_merged_named.pvals.gz"
ROADMAP_META="${ANNOTDIR}/${GENOME}/eid_to_mnemonic.txt"
PBC="${WORKDIR}/$SAMPLE.pbc.qc"


GENOME='hg19' # This is the only genome that currently works
SAMPLE="$Sample"

PBC="${WORKDIR}/$SAMPLE.pbc.qc"
FINAL_BAM="${Sample}/${Sample}.sorted.nodup.noM.bam"
FINAL_BED="${Sample}/${Sample}.nodup.tn5.tagAlign.gz"

F1="${Sample}.R1.trim.fastq.gz"
F2="${Sample}.R2.trim.fastq.gz"
ALIGNED_BAM="${Sample}/${Sample}.sorted.bam"
WORKDIR=$PWD
OUTDIR="qc"
OUTPREFIX=$Sample
INPREFIX=$Sample
GENOME='hg19' # This is the only genome that currently works

# Annotation files
X="/athena/elementolab/scratch/asd2007/Reference/hg19/"
DNASE_BED="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz"
BLACKLIST_BED="${ANNOTDIR}/${GENOME}/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
TSS_BED="${ANNOTDIR}/${GENOME}/hg19_RefSeq_stranded.bed.gz"
REF_FASTA="${ANNOTDIR}/${GENOME}/encodeHg19Male.fa"
PROM="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_prom_p2.bed.gz"
ENH="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_enh_p2.bed.gz"
REG2MAP="${ANNOTDIR}/${GENOME}/dnase_avgs_reg2map_p10_merged_named.pvals.gz"
ROADMAP_META="${ANNOTDIR}/${GENOME}/eid_to_mnemonic.txt"


#cp "${Sample}/${Sample}.sorted.bam" "${Sample}/${Sample}.sort.bam"


#python /home/asd2007/ATACseq/run_ataqc.athena.py --workdir $PWD/${Sample} \


python /home/asd2007/ATACseq/run_ataqc.athena.py --workdir $PWD/${Sample} \
    --outdir qc \
    --outprefix ${Sample} \
    --genome hg19 \
    --ref "${REFDIR}/hg19/encodeHg19Male/encodeHg19Male.fa" \
    --tss "${REFDIR}/hg19/hg19_RefSeq_stranded.bed.gz" \
    --dnase "${REFDIR}/hg19/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz" \
    --blacklist "${REFDIR}/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz" \
    --prom "${REFDIR}/hg19/reg2map_honeybadger2_dnase_prom_p2.bed.gz" \
    --enh "${REFDIR}/hg19/reg2map_honeybadger2_dnase_enh_p2.bed.gz" \
    --reg2map ${REG2MAP} \
    --meta ${ROADMAP_META} \
    --fastq1 ${F1} \
    --fastq2 ${F2} \
    --alignedbam ${ALIGNED_BAM} \
    --alignmentlog "${Sample}/${Sample}.align.log" \
    --coordsortbam "${Sample}/${Sample}.sorted.bam" \
    --duplog "${Sample}/${Sample}.dup.qc" \
    --pbc "${Sample}/${Sample}.pbc.qc" \
    --finalbam "${FINAL_BAM}" \
    --finalbed "${FINAL_BED}" \
    --bigwig "$Sample/$Sample.smooth150.center.extend.fpkm.max150.bw" \
    --peaks "${Sample}/${Sample}.nodup.tn5.pval0.1.500k.narrowPeak.gz" \
    --naive_overlap_peaks "${Sample}/pseudo_reps/${Sample}.nodup.tn5.pooled.pf.pval0.1.500K.naive_overlap.narrowPeak.gz" \
    --idr_peaks "${Sample}/IDR/${Sample}.IDR.txt.IDR0.1.filt.narrowPeak.gz"

