#!/bin/bash -l
#$ -N ChipSeqAlign
#$ -j y
#$ -m a
#$ -cwd
#$ -l zenodotus=true
#$ -l os=rhel6.3
#$ -M ashley.doane@gmail.com
#$ -l h_rt=32:00:00
#$ -pe smp 4-8
#$ -l h_vmem=12G
#$ -R y
#$ -o /home/asd2007/joblogs

BT2=1
BWA=0
INPUT_FASTQ=1
path=$1 #path to all the Sample folders
RSA=0
RSAINPUT=0
#file=`awk 'NR==n' n=$SGE_TASK_ID INPUTS`
gtf_path=$2 #Number indicating the reference gtf [1 for Human 2 for Mouse 3 for other]
#CONT="/home/asd2007/melnick_bcell_scratch/asd2007/Projects/barbieri/PCa/Zhao/sample_SRR2033042_LNCaP_input/sample_SRR2033042_LNCaP_input/sample_SRR2033042_LNCaP_input.sorted.nodup.bam"
#CONT=$3
#CONT="/home/asd2007/melnick_bcell_scratch/asd2007/ChipSeq/GCB_PolII/Sample_GCB_INPUT_C/Sample_GCB_INPUT_C/Sample_GCB_INPUT_C.tagAlign.gz"



#CONT="/home/asd2007/melnick_bcell_scratch/asd2007/ChipSeq/GCB_PolII/Sample_GCB_INPUT_C/Sample_GCB_INPUT_C_picard.bam"

# Uses job array for each Sample in the folder
file=$(ls ${path} | tail -n +${SGE_TASK_ID}| head -1)

#change directory to node directory
cd ${TMPDIR}

echo "File : $file Read from $path\n"

#Obtain the name of the Sample
Sample=$(basename "$file")
Sample=${Sample%%.*}

#Add file checker if fastq or gz

#copy the Samples with .gz extension to the tmp folder
#rsync -r -v -a -z $path/$file/*.fastq ./

#rsync -r -v -a -z $path/$file/*.fastq ./

rsync  -a $path/$file/*.sra ./
rsync  -a  $path/$file/*.gz ./

rsync -r -a   $path/$file/ ./


##for SRA

echo "making fastq fromm SRA using fastq-dump.."

#fastq-dump *.sra
mkdir -p INPUT
mv *INPUT.* INPUT/


if [ $RSAINPUT == 1 ]
then
    fastq-dump --gzip --skip-technical  --readids --dumpbase --split-files --clip *INPUT*.sra
fi
#rm *INPUT*.sra

myd=$(ls -lrth $TMPDIR)
echo "Dir contains $myd"
#rsync -r -v -a -z $path/$file/ ./
#rsync -r -v -a -z $path/$file/ ./

echo "Processing  $Sample ..."

#Figuring out the Reference genome

if [ $gtf_path == 1 ]
then
	#gtf="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/mm10_UCSC_ref.gtf"
	#REF="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex"
	REF="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome"
  REFbt2="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index"
   # REFbt="/home/asd20i07/dat02/asd2007/Reference/Homo_sapiens/UCSC/mm10/Sequence/BowtieIndex/genome"
  #  BLACK="/home/asd2007/dat02/asd2007/Reference/encodeBlack.bed"
  BLACK="/home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/Anshul_Hg19UltraHighSignalArtifactRegions.bed"
  fetchChromSizes hg19 > hg19.chrom.sizes
    chrsz="hg19.chrom.sizes"
    "/home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19.genome.chrom.sizes"
    RG="hg19"
elif [ $gtf_path == 2 ]
then
	gtf="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Mus_UCSC_ref.gtf"
	REF="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex"
    REFbt2="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index"
    BLACK="/home/asd2007/dat02/asd2007/Reference/mm10-blacklist.bed"
    RG="mm10"
    REFGen="/home/asd2007/melnick_bcell_scratch/asd2007/bin/bcbio/genomes/Mmusculus/mm10/seq/"
    chrsz="/home/asd2007/melnick_bcell_scratch/asd2007/Reference/mm10.genome.chrom.sizes"
else
	gtf="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/mm10_UCSC_ref.gtf"
	REF="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Homo_sapiens/UCSC/mm10/Sequence/BWAIndex/genome.fa"
fi

#rsync -a ${REF} ${TMPDIR}/

rsync -a ${REFbt2} ${TMPDIR}/

mkdir ${TMPDIR}/${Sample}

echo "ls of pwd is"
ls -lrth

ls -lrth $TMPDIR

#rsync -avP /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Reference/bwaidx/genome.fa* ref

#rsync -avP /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome
 #rsync -avP /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/ ref
#REF="$TMPDIR/ref/genome.fa"

mkdir R1



#mv *_1.fastq R1/
#mv *_2.fastq R2/



#count=$(ls -l R1/ |wc -l)
echo "----------bt2 aligning-------------"
#Use bwa mem to align each pair in data set


ReadGroup="$( echo '@RG'"\tID:${Sample}\tLB:${Sample}\tPL:illumina\tSM:${Sample}\tPU:${Sample}")"
if [ $INPUT_FASTQ == 1 ]
then
    mkdir -p $TMPDIR/INPUT
    cat ${TMPDIR}/*INPUT*.fastq.gz > ${TMPDIR}/INPUT/${Sample}.INPUT.fastq.gz
    rm ${TMPDIR}/*INPUT*.fastq.gz
    bwa mem -t ${NSLOTS} ${REF}*.fa $TMPDIR/INPUT/${Sample}.INPUT.fastq.gz >  $TMPDIR/${Sample}.INPUT.sam 
   # bowtie2  --very-sensitive-local -p ${NSLOTS} -x $TMPDIR/Bowtie2Index/genome -U $TMPDIR/INPUT/${Sample}.INPUT.fastq.gz -S $TMPDIR/${Sample}.INPUT.sam
    samtools view -bS $TMPDIR/${Sample}.INPUT.sam | samtools sort  - -o $TMPDIR/${Sample}/${Sample}.INPUT.sorted.bam
    rm *.sam
    samtools index $TMPDIR/${Sample}/${Sample}.INPUT.sorted.bam
    mkdir -p ${TMPDIR}/tmp
    CONT=$TMPDIR/${Sample}/${Sample}.INPUT.sorted.bam
    java -Xmx8g -jar /home/ole2001/PROGRAMS/SOFT/picard-tools-1.71/MarkDuplicates.jar REMOVE_DUPLICATES=true INPUT=$CONT METRICS_FILE=$TMPDIR/${Sample}/Input.dedup.txt OUTPUT=$TMPDIR/${Sample}/${Sample}.INPUT.sorted.nodup.bam VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=True
    CONT="${TMPDIR}/${Sample}/${Sample}.INPUT.sorted.nodup.bam"
    samtools index $CONT
fi


#CONT="${TMPDIR}/${Sample}/${Sample}.INPUT.sorted.nodup.bam"

###########
echo "----- FINISHED PROCESSING INPUT------"


if [ $RSA == 1 ]
then
    fastq-dump --gzip --skip-technical  --readids --dumpbase --split-files --clip *.sra
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
        bwa mem -t ${NSLOTS} ${REF}*.fa $TMPDIR/R1/${Sample}.fastq.gz | \
        samtools view -bS - > $TMPDIR/${Sample}/${Sample}.bam
    sambamba sort --memory-limit 35GB \
             --nthreads ${NSLOTS} --tmpdir ${TMPDIR} --out $TMPDIR/${Sample}/${Sample}.sorted.bam $TMPDIR/${Sample}/${Sample}.bam
    samtools index $TMPDIR/${Sample}/${Sample}.sorted.bam
fi
#bwa mem -t ${NSLOTS} ${REF}*.fa $TMPDIR/R1/${Sample}.fq > ${TMPDIR}/${Sample}/${Sample}.sam




echo "----------------------Remove duplicates--------------------"
#Remove duplicates using picard
#need to fix/rebuild picard

mkdir -p ${TMPDIR}/tmp
picard MarkDuplicates  INPUT=$TMPDIR/${Sample}/${Sample}.sorted.bam METRICS_FILE=$TMPDIR/${Sample}/dedup.txt REMOVE_DUPLICATES=true OUTPUT=$TMPDIR/${Sample}/${Sample}.sorted.nodup.bam VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=True



samtools index $TMPDIR/${Sample}/${Sample}.sorted.nodup.bam



#java -Xmx4g -jar /home/akv3001/Programs/picard/dist/picard.jar MarkDuplicates INPUT=$TMPDIR/${Sample}/${Sample}.sorted.bam METRICS_FILE=$TMPDIR/${Sample}/dedup.txt \
#    OUTPUT=$TMPDIR/${Sample}/${Sample}.sorted.noDUPS.bam VALIDATION_STRINGENCY=LENIENT \
#    ASSUME_SORTED=true REMOVE_DUPLICATES=true TMP_DIR=${TMPDIR}/tmp


#rsync -r -v $TMPDIR/${Sample} $path/${Sample}



echo "--------------------------Removing encode black listed intervals---------------------------"

bedtools subtract -A -a $TMPDIR/${Sample}/${Sample}.sorted.nodup.bam -b $BLACK > $TMPDIR/${Sample}/${Sample}.sorted.nodup.black.bam

samtools index $TMPDIR/${Sample}/${Sample}.sorted.nodup.black.bam



#/home/ole2001/PROGRAMS/SOFT/bedtools2/bin/bamToBed -i $TMPDIR/${Sample}/${Sample}.sorted.nodup.black.bam  > $TMPDIR/${Sample}/${Sample}.sorted.nodup.black.bed


rsync -a -r $TMPDIR/${Sample} $path/${Sample}


echo "------------------------------------------Call Peaks with MACS2--------------------------------------------"

#CONT="/home/asd2007/melnick_bcell_scratch/asd2007/ChipSeq/MD901shCrebpPooledInput/Sample_md901_pooled_input/Sample_md901_pooled_input/Sample_md901_pooled_input.sorted.bam"
#CONT="/home/asd2007/melnick_bcell_scratch/asd2007/ChipSeq/GCB_INPUT/GCB_INPUT.bt2.sorted.rmdup.bam"
#CONT="/home/asd2007/dat02/asd2007/Projects/barbieri/chipSeq/Sample_I_1/Sample_I_1/Sample_I_1.hg19.sorted.nodup.bam"
#rsync -avP /home/asd2007/melnick_bcell_scratch/asd2007/ChipSeq/MD901shCrebpPooledInput/Sample_md901_pooled_input/Sample_md901_pooled_input/Sample_md901_pooled_input.sorted.bam $TMPDIR/
#rsync -avP /home/asd2007/dat02/asd2007/Projects/barbieri/chipSeq/Sample_I_1/Sample_I_1/Sample_I_1.sorted.nodup.black.bam $TMPDIR/${Sample}/


#CONT=${Sample}.INPUT.sorted.bam

#rsync -avP ${CONT} $TMPDIR/

#CONT="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/ChIPSeq/LY7_wendy/LY7_669/Sample_Ly7_Input_669/Ly7_INPUT_669.INPUT.sorted.nodup.bam"
samtools view -b ${CONT} | bamToBed -i stdin | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c > $TMPDIR/Input.tagAlign.gz

samtools view -b $TMPDIR/${Sample}/${Sample}.sorted.nodup.bam | bamToBed -i stdin | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c > $TMPDIR/${Sample}/${Sample}.tagAlign.gz
#rsync -avP /home/asd2007/dat02/asd2007/Projects/barbieri/chipSeq/Sample_I_1/Sample_I_1/Sample_I_1.sorted.nodup.black.bam $TMPDIR/
cp $TMPDIR/Input.tagAlign.gz $TMPDIR/${Sample}/${Sample}.Input.tagAlign.gz

#macs2 callpeak -t $TMPDIR/${Sample}/${Sample}.tagAlign.gz -c $TMPDIR/Input.tagAlign.gz -f BED -n $TMPDIR/${Sample}/${Sample}.narrow -g hs -p 1e-3 --keep-dup all --call-summits --to-large --bdg



macs2 callpeak -t $TMPDIR/${Sample}/${Sample}.tagAlign.gz -c $TMPDIR/Input.tagAlign.gz -f BED -n $TMPDIR/${Sample}/${Sample}.narrow -g hs -p 1e-3 --keep-dup all --call-summits --to-large -B --SPMR



#macs2 callpeak -t $TMPDIR/${Sample}/${Sample}.tagAlign.gz -c $TMPDIR/Input.tagAlign.gz -f BED -n $TMPDIR/${Sample}/${Sample}.broad -g hs --keep-dup all --broad --broad-cutoff 0.1 --to-large



rsync -r -a  $TMPDIR/${Sample} $path/${Sample}/



fetchChromSizes hg19 > hg19.chrom.sizes

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

rsync -a ${prefix}_fc.bw /home/asd2007/melnick_bcell_scratch/asd2007/DATASETS/chipSeqSCORE/

#slopBed -i ${prefix}_FE.bdg -g $chrsz -b 0 | bedClip stdin $chrsz ${prefix}.FE.bedGraph

#bedGraphToBigWig ${p1}.FE.bedGraph ${chrsz} ${prefix}.FE.bigWig

#sort -k 1,1 -k 2,2n ${prefix}_fc.bdg > ${prefix}_fc_srt.bdg
#bedGraphToBigWig ${prefix}_fc_srt.bdg $chrsz ${prefix}_fc.bigwig
#rm -f ${prefix}_fc_srt.bdg ${prefix}_fc.bdg

#rm ${p1}_FE.bdg ${p1}_treat_pileup.bdg ${p1}_control_lambda.bdg

# Create bigwig file
sort -k8nr,8nr ${prefix}_peaks.narrowPeak | gzip -c > ${prefix}.peaks.bed.gz
#rm ${p1}_peaks.encodePeak
zcat ${prefix}.peaks.bed.gz | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$8}' > ${prefix}.4col.peaks.bed


rsync -r -a $TMPDIR/${Sample} $path/${Sample}/



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




rsync -a -r $TMPDIR/${Sample} $path/${Sample}/
#
#
#


















#
#
