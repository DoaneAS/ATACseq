#!/bin/bash -l
#$ -N ChipSeqAlign
#$ -j y
#$ -m a
#$ -cwd
#$ -l "athena=true"
#$ -M ashley.doane@gmail.com
#$ -l h_rt=42:00:00
#$ -pe smp 4
#$ -l h_vmem=5G
#$ -R y
#$ -o /home/asd2007/joblogs

## This script runs an SGE array job for chipseq processing and analysis. It takes as input a folder containing a set of folders, one for each sample.  Each sample folder contains a target chipseq fastq.gz and an input fastq.gz that must contain the string "INPUT".


BT2=0
BWA=1
INPUT_FASTQ=1
FOLDERPATH=$1 #FOLDERPATH to all the Sample folders
RSA=0
RSAINPUT=0
#FILEPATH=`awk 'NR==n' n=$SGE_TASK_ID INPUTS`
#gtf_FOLDERPATH=$2 #Number indicating the reference gtf [1 for Human 2 for Mouse 3 for other]
#CONT="/home/asd2007/melnick_bcell_scratch/asd2007/Projects/barbieri/PCa/Zhao/sample_SRR2033042_LNCaP_input/sample_SRR2033042_LNCaP_input/sample_SRR2033042_LNCaP_input.sorted.nodup.bam"
#CONT=$3
#CONT="/home/asd2007/melnick_bcell_scratch/asd2007/ChipSeq/GCB_PolII/Sample_GCB_INPUT_C/Sample_GCB_INPUT_C/Sample_GCB_INPUT_C.tagAlign.gz"



while [[ $# -gt 1 ]]
do
    key="$1"

    case $key in
        -f|--folder)
            FOLDERPATH="$2"
            shift # past argument
            ;;
        -g|--genome)
            GENOME="$2"
            shift # past argument
            ;;
        -t|--trim)
            TRIM=YES
            shift # past argument
            ;;
        --align)
            BWA=1
            ;;
        *)
            # unknown option
            ;;
    esac
    shift # past argument or value
done
echo FOLDER PATH   = "${FOLDERPATH}"
echo GENOME     = "${GENOME}"
echo TRIM READS    = "${TRIM}"
echo BWA ALN    = "${BWA}"
if [[ -n $1 ]]; then
    echo "Last line of FILEPATH specified as non-opt/last argument:"
    tail -1 $1
fi


BT2LAN=0

#set -e

# detect kernel then source spack FILEPATH if it exists

#if [ -f /pbtech_mounts/softlib001/apps/EL6/spack/share/spack/setup-env.sh ] ; then
#    . /pbtech_mounts/softlib001/apps/EL6/spack/share/spack/setup-env.sh
#fi


spack load bwa
spack load bowtie2
source /softlib/apps/EL6/R-v3.3.0/env



#if [ -f /home/asd2007/ATACseq/set-env.sh ] ; then
#     . /home/asd2007/ATACseq/set-env.sh

   


#FOLDERPATH=$1 #FOLDERPATH to all the Samples
#GENOME=$2
#GENOME="hg38" #Number indicating the reference gtf [1 for Human 2 for Mouse 3 for other]

# Uses job array for each Sample in the folder
FILEPATH=$(ls ${FOLDERPATH} | tail -n +${SGE_TASK_ID}| head -1)

#change directory to node directory
cd $TMPDIR

echo "File : $FILEPATH Read from $FOLDERPATH\n"

#Obtain the name of the Sample
Sample=$(basename "$FILEPATH")
Sample=${Sample%%.*}


rsync -r -v -a  $FOLDERPATH/$FILEPATH/*.gz ./
#rsync -r -v -a -z $FOLDERPATH/$FILEPATH/*.sra ./
rsync -r -v -a $FOLDERPATH/$FILEPATH/* ./
#SRA=$(ls *.sra)

mkdir -p ${Sample}

echo "Processing  $Sample ..."

#Figuring out the reference genome


if [[ $ATHENA == 1 ]] ; then
    REFDIR="/athena/elementolab/scratch/asd2007/reference"
    ANNOTDIR="/athena/elementolab/scratch/asd2007/reference"
    PICARDCONF="/home/asd2007/Scripts/picardmetrics.athena.conf"
    export PATH="/home/asd2007/anaconda2/bin:$PATH"
else
    REFDIR="/zenodotus/dat01/melnick_bcell_scratch/asd2007/reference"
    ANNOTDIR="/zenodotus/dat01/melnick_bcell_scratch/asd2007/reference"
    PICARDCONF="/home/asd2007/Scripts/picardmetrics.conf"
fi

# get genome args
if [[ $GENOME == "hg19" ]] ; then
    REF="/athena/elementolab/scratch/asd2007/reference/hg19/bwa_index/male.hg19.fa"
    #REF="${REFDIR}/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex"
	  REFbt2="${REFDIR}/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index"
    BLACK="/athena/elementolab/scratch/asd2007/reference/hg19/wgEncodeDacMapabilityConsensusExcludable.bed"
    #BLACK="${REFDIR}/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
    #fetchChromSizes hg19 > hg19.chrom.sizes
    chrsz="/athena/elementolab/scratch/asd2007/reference/hg19.chrom.sizes"
    cp $chrsz $PWD/hg19.chrom.sizes
    RG="hg19"
    SPEC="hs"
    REFGen="/athena/elementolab/scratch/asd2007/reference/hg19/seq/male.hg19.fa"
elif [[ $GENOME == "hg38" ]] ; then
    echo "genome is ${GENOME}"
    DNASE_BED="${ANNOTDIR}/${GENOME}/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz"
    BLACK="/athena/elementolab/scratch/asd2007/reference/hg38/hg38.blacklist.bed.gz"
    #PROM="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_prom_p2.bed.gz"
    #ENH="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_enh_p2.bed.gz"
    #REG2MAP="${ANNOTDIR}/${GENOME}/dnase_avgs_reg2map_p10_merged_named.pvals.gz"
    #ROADMAP_META="${ANNOTDIR}/${GENOME}/eid_to_mnemonic.txt"
	  REF="/athena/elementolab/scratch/asd2007/reference/hg38/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
	  REFbt2="/athena/elementolab/scratch/asd2007/reference/hg38/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
	  bwt2_idx="/athena/elementolab/scratch/asd2007/reference/hg38/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
    #fetchChromSizes hg19 > hg19.chrom.sizes
    RG="hg38"
    SPEC="hs"
    REFGen="/athena/elementolab/scratch/asd2007/reference/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
    chrsz=/athena/elementolab/scratch/asd2007/reference/hg38/hg38.chrom.sizes
    seq=/athena/elementolab/scratch/asd2007/reference/hg38/seq
    gensz=hs
    bwt2_idx=/athena/elementolab/scratch/asd2007/reference/hg38/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
    REF_FASTA=/athena/elementolab/scratch/asd2007/reference/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
    species_browser=hg38
    # data for ATAQC
    TSS_ENRICH=/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/hg38_gencode_tss_unique.bed.gz
    DNASE=/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.hg19_to_hg38.bed.gz
    PROM=/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/reg2map_honeybadger2_dnase_prom_p2.hg19_to_hg38.bed.gz
    ENH=/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/reg2map_honeybadger2_dnase_enh_p2.hg19_to_hg38.bed.gz
    REG2MAP=/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/hg38_dnase_avg_fseq_signal_formatted.txt.gz
    REG2MAP_BED=/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/hg38_celltype_compare_subsample.bed.gz
    ROADMAP_META=/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/hg38_dnase_avg_fseq_signal_metadata.txt
elif [[ $GENOME == "hg38_ENCODE" ]]; then
    DNASE_BED="${ANNOTDIR}/${GENOME}/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz"
    BLACK="/athena/elementolab/scratch/asd2007/reference/hg38/hg38.blacklist.bed.gz"
    #PROM="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_prom_p2.bed.gz"
    #ENH="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_enh_p2.bed.gz"
    #REG2MAP="${ANNOTDIR}/${GENOME}/dnase_avgs_reg2map_p10_merged_named.pvals.gz"
    #ROADMAP_META="${ANNOTDIR}/${GENOME}/eid_to_mnemonic.txt"
	  REF="/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
	  REFbt2="/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
    chrsz="/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/hg38_ENCODE.chrom.sizes"
    #fetchChromSizes hg19 > hg19.chrom.sizes
    RG="hg38"
    SPEC="hs"
    REFGen="/athena/elementolab/scratch/asd2007/reference/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
    chrsz=/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/hg38.chrom.sizes
    seq=/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/seq
    gensz=hs
    bwt2_idx=/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
    REF_FASTA=/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
    species_browser=hg38
    # data for ATAQC
    TSS_ENRICH=/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/ataqc/hg38_gencode_tss_unique.bed.gz
    DNASE=/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.hg19_to_hg38.bed.gz
    PROM=/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/ataqc/reg2map_honeybadger2_dnase_prom_p2.hg19_to_hg38.bed.gz
    ENH=/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/ataqc/reg2map_honeybadger2_dnase_enh_p2.hg19_to_hg38.bed.gz
    REG2MAP=/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/hg38_ENCODE_dnase_avg_fseq_signal_formatted.txt.gz
    REG2MAP_BED=/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/ataqc/hg38_celltype_compare_subsample.bed.gz
    ROADMAP_META=/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/ataqc/hg38_dnase_avg_fseq_signal_metadata.txt
elif [[ $GENOME == "mm10" ]]; then
    genome=mm10
	  gtf="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Mus_UCSC_ref.gtf"
	  REF="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex"
    REFbt2="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index"
    BLACK="/athena/elementolab/scratch/asd2007/reference/mm10-blacklist.bed"
    RG="mm10"
    SPEC="mm"
    REFGen="/athena/elementolab/scratch/asd2007/bin/bcbio/genomes/Mmusculus/mm10/seq/"
    #chrsz="/athena/elementolab/scratch/asd2007/reference/mm10.genome.chrom.sizes"
    fetchChromSizes mm10 > mm10.chrom.sizes
    chrsz= $PWD/mm10.chrom.sizes
    rsync -av /home/asd2007/Scripts/picardmetrics.Mouse.conf ./
else
    echo "genome is hg19"
	  REF="${REFDIR}/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex"
	  REFbt2="${REFDIR}/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index"
    BLACK="${REFDIR}/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
    #fetchChromSizes hg19 > hg19.chrom.sizes
    chrsz="/athena/elementolab/scratch/asd2007/reference/hg19.chrom.sizes"
    cp $chrsz $PWD/hg19.chrom.sizes
    RG="hg19"
    SPEC="hs"
    REFGen="/athena/elementolab/scratch/asd2007/local/share/bcbio/genomes/Hsapiens/hg19/seq/hg19.fa"
fi

# Uses job array for each Sample in the folder
##for SRA

echo "making fastq fromm SRA using fastq-dump.."

#fastq-dump *.sra


if [ $RSAINPUT == 1 ]
then
    fastq-dump --gzip --skip-technical  --readids --dumpbase --split-FILEPATHs --clip *INPUT*.sra
fi
#rm *INPUT*.sra

myd=$(ls -lrth $TMPDIR)
echo "Dir contains $myd"
#rsync -r -v -a -z $FOLDERPATH/$FILEPATH/ ./
#rsync -r -v -a -z $FOLDERPATH/$FILEPATH/ ./

echo "Processing  $Sample ..."
mkdir ${TMPDIR}/${Sample}

echo "ls of pwd is"
ls -lrth $TMPDIR
mkdir R1

#count=$(ls -l R1/ |wc -l)

ReadGroup="$( echo '@RG'"\tID:${Sample}\tLB:${Sample}\tPL:illumina\tSM:${Sample}\tPU:${Sample}")"
if [ $INPUT_FASTQ == 1 ]
then
    echo "----------bwa aligning INPUT-------------"
    #Use bwa mem to align each pair in data set
    mkdir -p $TMPDIR/INPUT
    cat ${TMPDIR}/*INPUT*.fastq.gz > ${TMPDIR}/INPUT/${Sample}.INPUT.fastq.gz
    #rm input from wd leaving only target
    rm ${TMPDIR}/*INPUT*.fastq.gz
    bwa mem -t ${NSLOTS} ${REF} $TMPDIR/INPUT/${Sample}.INPUT.fastq.gz | samtools view -bS - > $TMPDIR/${Sample}/${Sample}.INPUT.bam
    sambamba sort --memory-limit 35GB \
             --nthreads ${NSLOTS} --tmpdir ${TMPDIR} --out $TMPDIR/${Sample}/${Sample}.INPUT.sorted.bam $TMPDIR/${Sample}/${Sample}.INPUT.bam
    samtools index $TMPDIR/${Sample}/${Sample}.INPUT.sorted.bam
    mkdir -p ${TMPDIR}/tmp
    CONTA=$TMPDIR/${Sample}/${Sample}.INPUT.sorted.bam
    picard MarkDuplicates REMOVE_DUPLICATES=true INPUT=$CONTA METRICS_FILE=$TMPDIR/${Sample}/Input.dedup.txt OUTPUT=${TMPDIR}/${Sample}/${Sample}.INPUT.sorted.nodup.bam VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=True
    CONT=${TMPDIR}/${Sample}/${Sample}.INPUT.sorted.nodup.bam
    samtools index $CONT
fi


#CONT="${TMPDIR}/${Sample}/${Sample}.INPUT.sorted.nodup.bam"

###########
echo "----- FINISHED PROCESSING INPUT------"


if [ $RSA == 1 ]
then
    fastq-dump --gzip --skip-technical  --readids --dumpbase --split-FILEPATHs --clip *.sra
fi

cat $TMPDIR/*.fastq.gz > $TMPDIR/R1/${Sample}.fastq.gz


echo "----------bowtie2 aligning-------------"
#also works
#bowtie2 -p ${NSLOTS} --local -x ${TMPDIR}/bowtie2/hg19 -U $TMPDIR/R1/${Sample}.fq -S $TMPDIR/${Sample}/${Sample}.sam
#bowtie2  --very-sensitive-local -p ${NSLOTS} -x $TMPDIR/Bowtie2Index/genome -U $TMPDIR/R1/${Sample}.fastq.gz -S $TMPDIR/${Sample}.sam


if [ $BT2 == 1 ]
then
    bowtie2 --very-sensitive-local --threads ${NSLOTS} -x $TMPDIR/Bowtie2Index/genome \
        -U $TMPDIR/R1/${Sample}.fastq.gz  2> ${Sample}/${Sample}.align.log| \
        samtools view -bS - > $TMPDIR/${Sample}/${Sample}.bam
    sambamba sort --memory-limit 35GB \
             --nthreads ${NSLOTS} --tmpdir ${TMPDIR} --out $TMPDIR/${Sample}/${Sample}.sorted.bam $TMPDIR/${Sample}/${Sample}.bam
    samtools index $TMPDIR/${Sample}/${Sample}.sorted.bam
fi

if [ $BWA == 1 ]
then
        bwa mem -t ${NSLOTS} ${REF} $TMPDIR/R1/${Sample}.fastq.gz | \
        samtools view -bS - > $TMPDIR/${Sample}/${Sample}.bam
    sambamba sort --memory-limit 35GB \
             --nthreads ${NSLOTS} --tmpdir ${TMPDIR} --out $TMPDIR/${Sample}/${Sample}.sorted.bam $TMPDIR/${Sample}/${Sample}.bam
    samtools index $TMPDIR/${Sample}/${Sample}.sorted.bam
fi
#bwa mem -t ${NSLOTS} ${REF}*.fa $TMPDIR/R1/${Sample}.fq > ${TMPDIR}/${Sample}/${Sample}.sam




echo "----------------------Remove duplicates--------------------"

mkdir -p ${TMPDIR}/tmp
picard MarkDuplicates  INPUT=$TMPDIR/${Sample}/${Sample}.sorted.bam METRICS_FILE=$TMPDIR/${Sample}/dedup.txt REMOVE_DUPLICATES=true OUTPUT=$TMPDIR/${Sample}/${Sample}.sorted.nodup.bam VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=True



samtools index $TMPDIR/${Sample}/${Sample}.sorted.nodup.bam



#java -Xmx4g -jar /home/akv3001/Programs/picard/dist/picard.jar

#    OUTPUT=$TMPDIR/${Sample}/${Sample}.sorted.noDUPS.bam VALIDATION_STRINGENCY=LENIENT \
#    ASSUME_SORTED=true REMOVE_DUPLICATES=true TMP_DIR=${TMPDIR}/tmp


#rsync -r -v $TMPDIR/${Sample} $FOLDERPATH/${Sample}



echo "--------------------------Removing encode black listed intervals---------------------------"

rmBlack.sh $CONT $BLACK
CONT=$TMPDIR/${Sample}/${Sample}*INPUT*black.bam

rmBlack.sh $TMPDIR/${Sample}/${Sample}.sorted.nodup.bam $BLACK


#rsync -a -r $TMPDIR/${Sample} $FOLDERPATH/${Sample}


echo "------------------------------------------Call Peaks with MACS2--------------------------------------------"

#CONT="/home/asd2007/melnick_bcell_scratch/asd2007/ChipSeq/MD901shCrebpPooledInput/Sample_md901_pooled_input/Sample_md901_pooled_input/Sample_md901_pooled_input.sorted.bam"
#CONT="/home/asd2007/melnick_bcell_scratch/asd2007/ChipSeq/GCB_INPUT/GCB_INPUT.bt2.sorted.rmdup.bam"
#CONT="/home/asd2007/dat02/asd2007/Projects/barbieri/chipSeq/Sample_I_1/Sample_I_1/Sample_I_1.hg19.sorted.nodup.bam"
#rsync -avP /home/asd2007/melnick_bcell_scratch/asd2007/ChipSeq/MD901shCrebpPooledInput/Sample_md901_pooled_input/Sample_md901_pooled_input/Sample_md901_pooled_input.sorted.bam $TMPDIR/
#rsync -avP /home/asd2007/dat02/asd2007/Projects/barbieri/chipSeq/Sample_I_1/Sample_I_1/Sample_I_1.sorted.nodup.black.bam $TMPDIR/${Sample}/


#CONT=${Sample}.INPUT.sorted.bam

#rsync -avP ${CONT} $TMPDIR/

#CONT="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/ChIPSeq/LY7_wendy/LY7_669/Sample_Ly7_Input_669/Ly7_INPUT_669.INPUT.sorted.nodup.bam"
#CONT="/home/asd2007/zeno/Projects/DataSets/ChIPSeq/LY7_wendy/LY7_669/Sample_Ly7_Input_669/Ly7_INPUT_669.INPUT.sorted.nodup.bam"
#CONT="/athena/elementolab/scratch/asd2007/Projects/DataSets/ChIPSeq/LY7_wendy/LY7_343/Sample_Ly7_Input_343/"






samtools view -b ${CONT} | bamToBed -i stdin | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c > $TMPDIR/Input.tagAlign.gz

samtools view -b $TMPDIR/${Sample}/${Sample}.sorted.nodup.black.bam | bamToBed -i stdin | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c > $TMPDIR/${Sample}/${Sample}.tagAlign.gz

#rsync -avP /home/asd2007/dat02/asd2007/Projects/barbieri/chipSeq/Sample_I_1/Sample_I_1/Sample_I_1.sorted.nodup.black.bam $TMPDIR/

cp $TMPDIR/Input.tagAlign.gz $TMPDIR/${Sample}/${Sample}.Input.tagAlign.gz

#macs2 callpeak -t $TMPDIR/${Sample}/${Sample}.tagAlign.gz -c $TMPDIR/Input.tagAlign.gz -f BED -n $TMPDIR/${Sample}/${Sample}.narrow -g hs -p 1e-3 --keep-dup all --call-summits --to-large --bdg



macs2 callpeak -t $TMPDIR/${Sample}/${Sample}.tagAlign.gz -c $TMPDIR/Input.tagAlign.gz -f BED -n $TMPDIR/${Sample}/${Sample}.narrow -g hs -p 1e-3 --keep-dup all --call-summits --to-large -B --SPMR



#macs2 callpeak -t $TMPDIR/${Sample}/${Sample}.tagAlign.gz -c $TMPDIR/Input.tagAlign.gz -f BED -n $TMPDIR/${Sample}/${Sample}.broad -g hs --keep-dup all --broad --broad-cutoff 0.1 --to-large





cp /athena/elementolab/scratch/asd2007/Reference/hg19.chrom.sizes ./
#fetchChromSizes hg19 > hg19.chrom.sizes

prefix=$TMPDIR/${Sample}/${Sample}.narrow

#macs2 bdgcmp -t ${prefix}_treat_pileup.bdg -c ${prefix}_control_lambda.bdg -o ${prefix}_FE.bdg -m FE


macs2 bdgcmp -t ${prefix}_treat_pileup.bdg -c ${prefix}_control_lambda.bdg \
    --o-prefix ${prefix} -m FE

slopBed -i ${prefix}_FE.bdg -g hg19.chrom.sizes -b 0 | bedClip stdin hg19.chrom.sizes ${prefix}_fc.bedgraph
rm -f ${prefix}_FE.bdg

sort -k1,1 -k2,2n ${prefix}_fc.bedgraph > ${prefix}_fc.srt.bedgraph

bedGraphToBigWig  ${prefix}_fc.srt.bedgraph hg19.chrom.sizes ${prefix}_fc.bw

rm ${Sample}/*.bedgraph


rm ${prefix}*.bedgraph

#rsync -a ${prefix}_fc.bw /zenodotus/dat01/melnick_bcell_scratch/asd2007/DATASETS/chipSeqSCORE

#slopBed -i ${prefix}_FE.bdg -g $chrsz -b 0 | bedClip stdin $chrsz ${prefix}.FE.bedGraph

#bedGraphToBigWig ${p1}.FE.bedGraph ${chrsz} ${prefix}.FE.bigWig

#sort -k 1,1 -k 2,2n ${prefix}_fc.bdg > ${prefix}_fc_srt.bdg
#bedGraphToBigWig ${prefix}_fc_srt.bdg $chrsz ${prefix}_fc.bigwig
#rm -f ${prefix}_fc_srt.bdg ${prefix}_fc.bdg

#rm ${p1}_FE.bdg ${p1}_treat_pileup.bdg ${p1}_control_lambda.bdg

# Create bigwig FILEPATH
sort -k8nr,8nr ${prefix}_peaks.narrowPeak | gzip -c > ${prefix}.peaks.bed.gz
#rm ${p1}_peaks.encodePeak
zcat ${prefix}.peaks.bed.gz | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$8}' > ${prefix}.4col.peaks.bed





echo "---------------------------------- Signal Tracks ----------------------------"

#bamCoverage --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.bam --binSize 5 \
#    --outFileFormat bigwig --smoothLength 200  \
#    --normalizeTo1x 2150570000 \
#    --ignoreForNormalization chrX \
#    -o $TMPDIR/${Sample}/${Sample}.smooth.bw \
#    --centerReads --extendReads 200 --numberOfProcessors ${NSLOTS}


#bamCoverage -b $TMPDIR/${Sample}/${Sample}.sorted.nodup.black.bam -o $TMPDIR/${Sample}/${Sample}.sorted.nodup.black.smooth.bw --binSize 5 --outFileFormat bigwig --smoothLength 200 --normalizeUsingRPKM --centerReads --numberOfProcessors ${NSLOTS} --extendReads


bamCoverage --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.bam --binSize 5 \
    --outFileFormat bigwig --smoothLength 250  \
    --normalizeUsingRPKM \
    --ignoreForNormalization chrX \
    -o $TMPDIR/${Sample}/${Sample}.rpkm.smooth.bw \
    --centerReads --extendReads 125 --numberOfProcessors ${NSLOTS}



## cleanup

if [ -f "${Sample}.sorted.nodup.black.bam"] ; then
    rm *.sorted.bam
    rm *.sorted.nodup.bam
fi




rsync -a -r $TMPDIR/${Sample} $FOLDERPATH/${Sample}/
#
#
#


















#
#
