#!/bin/bash -l
#$ -N AtacSEQ
#$ -j y
#$ -m a
#$ -cwd
#$ -l zenodotus=true
#$ -l os=rhel6.3
#$ -M ashley.doane@gmail.com
#$ -l h_rt=40:50:00
#$ -pe smp 4-8
#$ -l h_vmem=12G
#$ -R y


path=$1 #path to all the Samples
gtf_path=$2 #Number indicating the reference gtf [1 for Human 2 for Mouse 3 for other]

echo "requirements: deeptools, python 2.7, macs2, picard, java, sambamba, bwa, ucsc-tools"

# Uses job array for each Sample in the folder
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

rsync -r -v -a -z $path/$file/*.gz ./

rsync -r -v -a -z $path/$file/*.sra ./
rsync -r -v -a -z $path/$file/ ./
#rsync -r -v -a -z $path/$file/ ./
#rsync -r -v -a -z $path/$file/ ./

#SRA=$(ls *.sra)


parallel -j ${NSLOTS} 'fastq-dump --split-files {}' ::: *.sra
#fastq-dump --split-files ${SRA}

#fastq-dump --split-files *.sra

echo "Processing  $Sample ..."

#Figuring out the Reference genome

if [ $gtf_path == 1 ]
then
	#gtf="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/mm10_UCSC_ref.gtf"
	REF="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex"
	REFbt2="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index"
   # REFbt="/home/asd20i07/dat02/asd2007/Reference/Homo_sapiens/UCSC/mm10/Sequence/BowtieIndex/genome"
    BLACK="/home/asd2007/dat02/asd2007/Reference/encodeBlack.bed"
    BLACK="/home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/Anshul_Hg19UltraHighSignalArtifactRegions.bed"
    chrsz="/home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19.genome.chrom.sizes"
    RG="hg19"
elif [ $gtf_path == 2 ]
then
	gtf="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Mus_UCSC_ref.gtf"
	REF="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex"
    REFbt2="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index"
    BLACK="/home/asd2007/dat02/asd2007/Reference/mm10-blacklist.bed"
    RG="mm10"
    chrsz="/home/asd2007/melnick_bcell_scratch/asd2007/Reference/mm10.genome.chrom.sizes"
else
	gtf="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/mm10_UCSC_ref.gtf"
	REF="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Homo_sapiens/UCSC/mm10/Sequence/BWAIndex/genome.fa"
fi

mkdir ref

#rsync -avP /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Reference/bwaidx/genome.fa* ref
#rsync -avP /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Homo_sapiens/UCSC/mm10/Sequence/BWAIndex/genome
 #rsync -avP /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Homo_sapiens/UCSC/mm10/Sequence/BWAIndex/ ref
#REF="$TMPDIR/ref/genome.fa"

rsync -avP ${REF} $TMPDIR/

rsync -avP ${REFbt2} $TMPDIR/

mkdir R1
mkdir R2
mkdir ${Sample}

gunzip *.gz


mv *_R1* R1/
mv *_R2* R2/

mv *_1* R1/
mv *_2* R2/


#mv *val_1* R1/
#mv *val_2* R2/
#mv *_1.fastq R1/
#mv *_2.fastq R2/



#count=$(ls -l R1/ |wc -l)
echo "----------bwa-mem aligning-------------"
#Use bwa mem to align each pair in data set

F1=$(ls R1)
F2=$(ls R2)

echo "aligning : $F1 , $F2 using bwa-mem.."

#gunzip -c $TMPDIR/R1/$F1 > $TMPDIR/R1/${Sample}.R1.fastq
#gunzip -c $TMPDIR/R2/$F2 > $TMPDIR/R2/${Sample}.R2.fastq
cat $TMPDIR/R1/*.*q > $TMPDIR/R1/${Sample}.R1.fastq
cat $TMPDIR/R2/*.*q > $TMPDIR/R2/${Sample}.R2.fastq


cutadapt -m 5 -e 0.20 -a CTGTCTCTTATA -A CTGTCTCTTATA -o $TMPDIR/R1/${Sample}.R1.fq -p $TMPDIR/R2/${Sample}.R2.fq $TMPDIR/R1/${Sample}.R1.fastq $TMPDIR/R2/${Sample}.R2.fastq


#bowtie2  -X 2000 -p ${NSLOTS} -x $TMPDIR/Bowtie2Index/genome -1 $TMPDIR/R1/${Sample}.R1.fq -2 $TMPDIR/R2/${Sample}.R2.fq -S $TMPDIR/${Sample}/${Sample}.bt2.sam


echo "aligning : $TMPDIR/R1/${Sample}.R1.fq ,  $TMPDIR/R2/${Sample}.R2.fq using bwa-mem.."

bwa mem -t ${NSLOTS} -M $TMPDIR/BWAIndex/genome.fa $TMPDIR/R1/${Sample}.R1.fq $TMPDIR/R2/${Sample}.R2.fq > $TMPDIR/${Sample}/${Sample}.sam




echo "----------Samtools Converting to Sorted bam-----------------"

#samtools view -bS $TMPDIR/${Sample}/${Sample}.sam| samtools sort  - $TMPDIR/${Sample}/${Sample}.sorted

samtools view -bS $TMPDIR/${Sample}/${Sample}.sam | samtools sort  - -o $TMPDIR/${Sample}/${Sample}.sorted.bam

#TMPDIR/${Sample}/${Sample}.bt2.sam

#samtools view -bS $TMPDIR/${Sample}/${Sample}.bt2.sam | samtools sort  - -o $TMPDIR/${Sample}/${Sample}.bt2.sorted.bam

echo "get insert Sizes of unfiltered alignments"

samtools view -f66 $TMPDIR/${Sample}/${Sample}.sorted.bam |cut -f 9|sed 's/^-//' > $TMPDIR/${Sample}/${Sample}.raw.InsertSizeMetrics.txt

#samtools view -f66 $TMPDIR/${Sample}/${Sample}.bt2.sorted.bam |cut -f 9|sed 's/^-//' > $TMPDIR/${Sample}/${Sample}.raw.bt2.InsertSizeMetrics.txt

#remove sam file after converted to sorted bam
rm $TMPDIR/${Sample}/${Sample}*.sam


echo "Create index"

samtools index  $TMPDIR/${Sample}/${Sample}.sorted.bam


RAW_BAM_FILE="$TMPDIR/${Sample}/${Sample}.sorted.bam"

#picard MarkDuplicates INPUT=${FILT_BAM_FILE} METRICS_FILE=$TMPDIR/${Sample}/dedup.txt \
#    OUTPUT=$TMPDIR/${Sample}.sorted.uniq.nodup.bam VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false
#RAW_BAM_FILE="$TMPDIR/${Sample}.sorted.uniq.nodup.bam"
# =============================
# Remove  unmapped, mate unmapped
# not primary alignment, reads failing platform
# Only keep properly paired reads
# ==================


mapq_thresh=30

FILT_BAM_PREFIX="${Sample}.filt.srt"
FILT_BAM_FILE="${FILT_BAM_PREFIX}.bam"
TMP_FILT_BAM_PREFIX="tmp.${FILT_BAM_PREFIX}.nmsrt"
TMP_FILT_BAM_FILE="${TMP_FILT_BAM_PREFIX}.bam"

samtools view -F 1804 -f 2 -q $mapq_thresh -u $RAW_BAM_FILE | \
    sambamba sort -n -t ${NSLOTS} /dev/stdin -o ${TMP_FILT_BAM_FILE}
samtools fixmate -r ${TMP_FILT_BAM_FILE} ${TMP_FILT_BAM_FILE}.fixmate.bam

samtools view -F 1804 -f 2 -u  ${TMP_FILT_BAM_FILE}.fixmate.bam | sambamba sort -t ${NSLOTS} /dev/stdin -o ${FILT_BAM_FILE}




echo "----------------------Remove duplicates--------------------"
#Remove duplicates using picard
#need to fix/rebuild picard
#java -Xmx8g -jar /home/ole2001/PROGRAMS/SOFT/picard-tools-1.71/MarkDuplicates.jar REMOVE_DUPLICATES=true INPUT=$TMPDIR/${Sample}/${Sample}.sorted.bam METRICS_FILE=$TMPDIR/${Sample}/dedup.txt OUTPUT=$TMPDIR/${Sample}/${Sample}.sorted.nodup.bam VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=True



java -Xmx8g -jar /home/akv3001/Programs/picard/dist/picard.jar MarkDuplicates INPUT=${FILT_BAM_FILE} METRICS_FILE=$TMPDIR/${Sample}/dedup.txt \
    OUTPUT=$TMPDIR/${Sample}.sorted.uniq.nodup.bam VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false



#java -Xmx4g -jar /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/picard-tools-1.119/MarkDuplicates.jar INPUT=$TMPDIR/${Sample}/${Sample}_sorted.bam METRICS_FILE=$TMPDIR/${Sample}/dedup.txt OUTPUT=$TMPDIR/${Sample}/${Sample}_noDuplicates.bam VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true TMP_DIR=/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Irene_2_Samples/tmp

#rsync -a -vP $TMPDIR/${Sample} $path/${Sample}


samtools view -F 1804 -f 2 -b ${TMPDIR}/${Sample}.sorted.uniq.nodup.bam > $TMPDIR/${Sample}/${Sample}.sorted.nodup.bam
sambamba index -t ${NSLOTS} $TMPDIR/${Sample}/${Sample}.sorted.nodup.bam
sambamba flagstat -t ${NSLOTS} $TMPDIR/${Sample}/${Sample}.sorted.nodup.bam > $TMPDIR/${Sample}/${Sample}.mapqc.txt


rsync -av /home/asd2007/Scripts/picardmetrics.conf ./

mkdir metrics

#picardmetrics run -o $TMPDIR/metrics $TMPDIR/${Sample}/${Sample}.sorted.bam


pbc_qc="$TMPDIR/${Sample}/${Sample}.library_complexity.qc.txt"


dupmark_bam=$TMPDIR/${Sample}.sorted.uniq.nodup.bam

## TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
sambamba sort -t ${NSLOTS} -n  $TMPDIR/${Sample}/${Sample}.sorted.nodup.bam -o $TMPDIR/${Sample}/${Sample}.nsorted.nodup.bam


bedtools bamtobed -bedpe -i  $TMPDIR/${Sample}/${Sample}.nsorted.nodup.bam | \
    awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | \
    grep -v 'chrM' | sort | uniq -c | \
    awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1_m2}' > $pbc_qc


echo "-----------------------------Remove chrm------------------------------"


samtools idxstats $TMPDIR/${Sample}/${Sample}.sorted.nodup.bam | cut -f 1 | grep -v chrM | xargs samtools view -b $TMPDIR/${Sample}/${Sample}.sorted.nodup.bam > $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.bam

samtools index $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.bam


#picardmetrics run -o $TMPDIR/metrics $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.bam

#mkdir /zenodotus/dat01/melnick_bcell_scratch/asd2007/atacData/QCmetrics/${Sample}
#rsync -r -v $TMPDIR/metrics /zenodotus/dat01/melnick_bcell_scratch/asd2007/atacData/QCmetrics/${Sample}/

echo "-----------------  additional QC ---------------------"

rsync -av /home/asd2007/Scripts/picardmetrics.conf ./

mkdir metrics

picardmetrics run -o $TMPDIR/metrics $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.bam


#picardmetrics run -o $TMPDIR/metrics $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.bam

mkdir /zenodotus/dat01/melnick_bcell_scratch/asd2007/atacData/QCmetrics/${Sample}
rsync -r -v $TMPDIR/metrics /zenodotus/dat01/melnick_bcell_scratch/asd2007/atacData/QCmetrics/${Sample}/

# histogram file
for w in 1000 500 200
do
    picard CollectInsertSizeMetrics I=$TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.bam O="$TMPDIR/${Sample}/${Sample}.window${w}.hist_data" H="$TMPDIR/${Sample}/${Sample}.window${w}.hist_graph.pdf" W=${w}
done

echo "--------------------------Removing encode black listed intervals---------------------------"



bedtools subtract -A -a $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.bam -b $BLACK > $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam

samtools index $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam

samtools flagstat  $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam > $TMPDIR/${Sample}/${Sample}.finalBam.mapStats.txt

rsync -avP $TMPDIR/${Sample} $path/${Sample}


echo "--------------------------Tn5 adjusted bedfile for MACS2 peak calling---------------------------"

samtools view -H $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam | grep chr | grep -v chrM | /home/ole2001/PERL_SCRIPTS/columns.pl 1 2 | sed 's/SN://' | sed 's/LN://' > chrom.sizes

#/home/ole2001/PROGRAMS/SOFT/bedtools2/bin/bamToBed -i $TMPDIR/${Sample}/${Sample}.mm10.sorted.nodup.noM.black.bam | \
#     sh /home/asd2007/bin/adjustBedTn5.sh  > $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bed
bamToBed -i $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam | \
    awk 'BEGIN{OFS="\t"} $6=="+" { $2=$2+4; $3=$3 ; $4="N" ; print $0 } $6=="-"{ $2=$2; $3=$3-5; $4="N" ; print $0 }' | gzip -c > $TMPDIR/${Sample}/${Sample}.tn5.tagAlign.gz


bedtools bamtobed -bedpe -mate1 -i  $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam | gzip -c > $TMPDIR/${Sample}/${Sample}.tn5.bedpe.tagAlign.gz

echo "----------------------------partition reads of <= 115= less than 123 after corrections ----------------------"
samtools view -h ${TMPDIR}/${Sample}/${Sample}.sorted.nodup.noM.bam | awk '/^@/ || $9 <= 123 && $9 >= -123' | samtools view -b - > ${TMPDIR}/${Sample}/${Sample}.sorted.nodup.noM.frag123.bam


echo "------------------------------------------Call Peaks with MACS2--------------------------------------------"


#adjustedBed="/home/ole2001/PROGRAMS/SOFT/bedtools2/bin/slopBed -i $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bed -g sizes -l 75 -r -75 -s"

#macs2 callpeak -t <(${adjustedBed}) -f BED -n $TMPDIR/${Sample}/${Sample}.broad -g 2.7e9 -p 1e-3 --nomodel --shift 75 -B --SPMR --broad --keep-dup all

#macs2 callpeak -t <(${adjustedBed}) -f BED -n $TMPDIR/${Sample}/${Sample}.narrow -g mm -p 1e-3 --nomodel --shift 75 -B --SPMR --keep-dup all --call-summits


macs2 callpeak -t $TMPDIR/${Sample}/${Sample}.tn5.tagAlign.gz -f BED -n $TMPDIR/${Sample}/${Sample}.narrow -g hs --nomodel --shift -75 --extsize 150 --keep-dup all --call-summits -p 1e-2

macs2 callpeak -t $TMPDIR/${Sample}/${Sample}.tn5.tagAlign.gz -f BED -n $TMPDIR/${Sample}/${Sample}.broad -g hs  --nomodel --shift -100 --extsize 200 --keep-dup all --broad --broad-cutoff 0.1

#macs2 callpeak -t $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam -f BAMPE -n $TMPDIR/${Sample}/${Sample}.narrow -g hs --nomodel --shift -75 --extsize 150 --keep-dup all --call-summits -p 1e-3

#macs2 callpeak -t $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam -f BAMPE -n $TMPDIR/${Sample}/${Sample}.broad -g hs  --nomodel --shift -100 --extsize 200 --keep-dup all --broad --broad-cutoff 0.1


#macs2 callpeak -t ${TMPDIR}/${Sample}/${Sample}.sorted.nodup.noM.frag123.bam -f BAMPE -n $TMPDIR/${Sample}/${Sample}.frag123.narrow -g hs --nomodel --shift -75 --extsize 150 --keep-dup all --call-summits -p 1e-3

## lenient threshold for IDR
#macs2 callpeak -t $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam -f BAMPE -n $TMPDIR/${Sample}/${Sample}.narrow.p0.01 -g hs --nomodel --shift -75 --extsize 150 --keep-dup all --call-summits -p 1e-1

#macs2 callpeak -t $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam -f BAMPE -n $TMPDIR/${Sample}/${Sample}.broad.p0.1 -g hs  --nomodel --shift -100 --extsize 200 --keep-dup all --broad --broad-cutoff 0.1
#rsync -avP /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.* ./
#rsync -avP /home/asd2007/dat02/asd2007/Reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.* ./
prefix="$TMPDIR/${Sample}/${Sample}.narrow"

macs2 bdgcmp -t "$prefix"_treat_pileup.bdg -c "$prefix"_control_lambda.bdg
    --o-prefix "$prefix" -m FE
slopBed -i "$prefix"_FE.bdg -g "$chrsz" -b 0 | bedClip stdin "$chrsz" "$prefix"_fc.bdg
rm -f "$prefix"_FE.bdg

sort -k1,1 -k2,2n "$prefix"_fc.bdg > "$prefix"_fc_srt.bdg
bedGraphToBigWig "$prefix"_fc_srt.bdg "$chrsz" "$prefix"_fc.bigwig
rm -f "$prefix"_fc_srt.bdg "$prefix"_fc.bdg


#echo "--- deeptools coverage -rpkm --"

#bamCoverage --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam --binSize 5 \
#    --outFileFormat bigwig --smoothLength 200  \
#    --normalizeUsingRPKM \
#    -o $TMPDIR/${Sample}/${Sample}.bamCov.rpkm.bw --centerReads --numberOfProcessors ${NSLOTS}

echo "----------- compute fragment midpoint coverage per bp ------------------"
#pyatac cov --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam --out $TMPDIR/${Sample}/${Sample}.pyatac.cov --cores ${NSLOTS}

#pyatac cov --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam --out $TMPDIR/${Sample}/${Sample}.pyatac.125bp.cov --cores ${NSLOTS} --upper 125

bamCoverage --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam --binSize 5 \
    --outFileFormat bigwig --smoothLength 150 \
    --normalizeUsingRPKM \
    -o $TMPDIR/${Sample}/${Sample}.smooth151.center.extend.fpkm.bw --centerReads --extendReads --numberOfProcessors ${NSLOTS}


bamCoverage --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam --binSize 20 \
    --outFileFormat bigwig --smoothLength 150 \
    --normalizeUsingRPKM \
    -o $TMPDIR/${Sample}/${Sample}.bin20.smooth150.fpkm.bw --numberOfProcessors ${NSLOTS}

rsync -avP $TMPDIR/${Sample} $path/${Sample}
rsync -r -a -v $TMPDIR/${Sample} $path/${Sample}

echo "------------------------------------------Call Nucleosomes with NucleoATAC------------------------------"
#
 bedtools slop -i $TMPDIR/${Sample}/${Sample}.broad_peaks.broadPeak -g chrom.sizes -b 1000 > $TMPDIR/${Sample}/${Sample}.slop1k.bed

 ## save to run on pooled samples
 ## high resolution 1bp bin insertion density
 #pyatac ins --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam --out $TMPDIR/${Sample}/${Sample}.ins.smooth --smooth 151 --cores ${NSLOTS}


#
#rsync -avP /home/asd2007/dat02/asd2007/Reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.* ./
#nucleoatac run --bed $TMPDIR/${Sample}/${Sample}.slop1k.bed --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam --fasta $TMPDIR/genome.fa --out  $TMPDIR/${Sample}/${Sample} \
#    --write_all --cores ${NSLOTS}
#
rsync -avP $TMPDIR/${Sample} $path/${Sample}
rsync -r -a -v $TMPDIR/${Sample} $path/${Sample}
#

## additional tracks


#igvtools toTDF -z 10 $TMPDIR/${Sample}/${Sample}.ins.smooth.ins.bedgraph.gz $TMPDIR/${Sample}/${Sample}.ins.smooth.ins.tdf $RG


#mv $TMPDIR/${Sample}/${Sample}.ins.smooth.ins.bedgraph.gz $TMPDIR/${Sample}/${Sample}.ins.smooth.ins.bdg.gz
#gunzip $TMPDIR/${Sample}/${Sample}.ins.smooth.ins.bdg.gz


#~/SHELL_SCRIPTS/bdg2bw $TMPDIR/${Sample}/${Sample}.ins.smooth.ins.bdg chrom.sizes


#mv $TMPDIR/${Sample}/${Sample}.nucleoatac_signal.smooth.bedgraph.gz $TMPDIR/${Sample}/${Sample}.nucleoatac_signal.smooth.bdg.gz

#gunzip $TMPDIR/${Sample}/${Sample}.nucleoatac_signal.smooth.bdg.gz
#~/SHELL_SCRIPTS/bdg2bw  $TMPDIR/${Sample}/${Sample}.nucleoatac_signal.smooth.bdg chrom.sizes


#rm $TMPDIR/${Sample}/${Sample}*.bdg

###
#rsync -avP $TMPDIR/${Sample} $path/${Sample}
#rsync -r -a -v $TMPDIR/${Sample} $path/${Sample}
