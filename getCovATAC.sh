#!/bin/bash -l
#$ -N pyatacINS
#$ -j y
#$ -m a
#$ -cwd
#$ -M ashley.doane@gmail.com
#$ -l athena=true
#$ -l h_rt=13:50:00
#$ -pe smp 6-8
#$ -l h_vmem=5G
#$ -R y
#$ -o /home/asd2007/joblogs

path=$1

file=$(ls ${path} | tail -n +${SGE_TASK_ID}| head -1)


#change directory to node directory
cd $TMPDIR

echo "File : $file Read from $path\n"

#Obtain the name of the Sample
Sample=$(basename "$file")
Sample=${Sample%%.*}




#Add file checker if fastq or gz

#copy the Samples with .gz extension to the tmp folder
#rsync -r -v -a -z $path/$file/*.fastq ./

#rsync -r -v -a -z $path/$file/*.fastq ./

rsync  -v -a  $path/$file/${Sample}/*.sorted.nodup.noM.bam ./

rmBlack.sh ${Sample}.sorted.nodup.noM.bam 

#rsync -r -v -a -z $path/$file ./

#rsync -r -v -a -z $file/ ./

#mv *.black.bam ${Sample}.bam

#samtools index ${Sample}.bam

#atlas=$2 #bed file with windows in which to to count inserts

#Atlas=$(basename "$atlas")

#file=$1
#change directory to node directory

#Add file checker if fastq or gz

#copy the samples with .gz extension to the tmp folder
#rsync -r -v -a -z --exclude 'Summary' $path/$file/*.gz ./

#rsync -r -v $file ./
#mkdir /home/asd2007/melnick_bcell_scratch/asd2007/atacData/pyatacIns/${Sample}

echo "Processing  $Sample ..."

#Figuring out the Reference genome 

#atlas="/home/asd2007/melnick_bcell_scratch/asd2007/atlas.beds/atlas.SE.bed"

#rsync -avP ${atlas} $TMPDIR/


echo "---------- indexing bam file  -------------"
samtools index ${TMPDIR}/${Sample}.sorted.nodup.noM.black.bam
echo "----------counting Tn5 insertions -------------"
#Use bwa mem to align each pair in data set

#pyatac counts --bam ${Sample}.bam --bed ${atlas} --out ${Sample}.${Atlas}.ins.115.counts.txt --upper 125

#pyatac counts --bam ${Sample}.bam --bed ${atlas} --out ${Sample}.${Atlas}.ins.counts.txt
#pyatac counts --bam ${Sample}.bam --bed ${atlas} --out ${Sample}.${Atlas}.ins.500.counts.txt --upper 1000

#zcat $TMPDIR/*.gz > ${Sample}.fastq


## pyatac cov setup 
#samtools view -H ${Sample}.bam | grep chr | grep -v chrM | /home/ole2001/PERL_SCRIPTS/columns.pl 1 2 | sed 's/SN://' | sed 's/LN://' > chrom.sizes

#pyatac cov --bam ${Sample}.bam --out ${Sample} --cores ${NSLOTS}

#mv ${Sample}.cov.bedgraph.gz ${Sample}.cov.bdg.gz

#gunzip ${Sample}.cov.bdg.gz

#~/SHELL_SCRIPTS/bdg2bw ${Sample}.cov.bdg chrom.sizes
##

#rsync -r -v $TMPDIR/${Sample}*.bw /zenodotus/dat01/melnick_bcell_scratch/asd2007/COVERAGE/ATAC/coverage/test/


#BED="/athena/elementolab/scratch/asd2007/Projects/DataSets/atacData/atacLY7/WTidr.IDR0.1.filt.narrowPeak"
#BED="/athena/elementolab/scratch/asd2007/Projects/HOMER/Homer.input/LY7idr.atlas.peakSummits300.bed"

BED="/athena/elementolab/scratch/asd2007/Projects/HOMER/Homer.input/LY7idr.atlas.peakSummits300reduced.bed"


pyatac counts --bam ${TMPDIR}/${Sample}.sorted.nodup.noM.black.bam --bed ${BED} --out ${Sample}.ins
#bedtools multicov -split -bams $BAMS -bed $BED

rsync -a -v $TMPDIR/${Sample}* /athena/elementolab/scratch/asd2007/Projects/DataSets/atacData/atacLY7/bams/

LIBS=$(zcat ${Sample}.ins.counts.txt.gz  | awk '{sum +=$1} END  {printf sum}')

let KM=20000000 #1 million * bin size in kb
let LIBZ=$LIBS/$KM



bamCoverage --bam ${TMPDIR}/${Sample}.sorted.nodup.noM.black.bam --binSize 20 \
    --outFileFormat bigwig --smoothLength 150  \
    --scaleFactor $LIBZ \
    -o ${TMPDIR}/$Sample.readsInPeaks.bin20.centered.smooth.150.max150f.bw --numberOfProcessors ${NSLOTS} \
    --maxFragmentLength 150 \
    --centerReads --extendReads



#bamCoverage --bam ${TMPDIR}/${Sample}.bam --binSize 20 \
#    --outFileFormat bigwig --smoothLength 150  \
#    --scaleFactor $LIBZ \
#    -o ${TMPDIR}/$Sample.readsInPeaks.bin20.centered.smooth.150.max500f.bw --numberOfProcessors ${NSLOTS} \
#    --maxFragmentLength 500 \
#    --centerReads --extendReads


let KM=5000000 #1 million * bin size in kb
let LIBZ=$LIBS/$KM

bamCoverage --bam ${TMPDIR}/${Sample}.sorted.nodup.noM.black.bam --binSize 5 \
    --outFileFormat bigwig --smoothLength 150  \
    --scaleFactor $LIBZ \
    -o ${TMPDIR}/$Sample.readsInPeaks.bin5.centered.smooth.150.bw --numberOfProcessors ${NSLOTS} \
    --maxFragmentLength 150 \
    --centerReads --extendReads


#bamCoverage --bam ${TMPDIR}/${Sample}.bam --binSize 5 \
#    --outFileFormat bigwig --smoothLength 150  \
#    --scaleFactor $LIBZ \
#    -o ${TMPDIR}/$Sample.readsInPeaks.bin5.centered.smooth.150.max500.bw --numberOfProcessors ${NSLOTS} \
#    --maxFragmentLength 500 \
#    --centerReads --extendReads
    #pyatac ins --bam $bam --bed $bed --

#path="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Jon_bwa_mm10_output"

rsync -a -v $TMPDIR/${Sample}* /athena/elementolab/scratch/asd2007/Projects/DataSets/atacData/atacLY7/bams/  







# --centerReads --extendReads 200 \#
#--numberOfProcessors ${NSLOTS}
#pyatac ins --bam $bam --bed $bed --

#path="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Jon_bwa_mm10_output"

#rsync -r -v $TMPDIR/${Sample}*.bw /zenodotus/dat01/melnick_bcell_scratch/asd2007/COVERAGE/ATAC/test

