#!/bin/bash

# read from command line the unfiltered and unsortde bam file

p1=$1


#BLACK="/home/asd2007/melnick_bcell_scratch/asd2007/Reference/encodeBlack.bed"
# help
if [ -z "$p1"  ]
then
    echo "This will sort bam, remove mitochondria/duplicates and histogram insert size"
    echo "$(basename $0) <bamFile>"
    echo "<samFile>: Required bam input file"
    exit
fi

#SUB=$2
#if [ $SUB == 1 ]#
#then
    # convert to bam and sort
 #   echo "Sorting..."
   # out1prefix=$(echo ${p1} | sed 's/\.bam$//')
   # out1="${out1prefix}.sorted.bam"
  #  echo ${out1}
  #  samtools view -s 0.01 -S -u $p1 | samtools sort - > ${out1}
    #index
 #   samtools index $out1
    # samtools sort  $TMPDIR/${Sample}.bam -o $TMPDIR/${Sample}/${Sample}.sorted.bam
    #cp $TMPDIR/${Sample}/${Sample}.sorted.bam  $TMPDIR/${Sample}/${Sample}.bt2.sorted.bam``
#else

    echo "Sorting..."
    out1prefix=$(echo $p1 | sed 's/\.bam$//')
    out1="${out1prefix}.sorted.bam"
    echo ${out1}
    #samtools view -S -u $p1 | samtools sort - > ${out1}
    #index
#    sambamba sort -t ${NSLOTS} $p1 > ${out1}

    sambamba sort --memory-limit 30GB \
             --nthreads ${NSLOTS} --tmpdir ${TMPDIR} --out ${out1} $p1
    samtools index $out1
    # echo "aligning : $TMPDIR/${Sample}.R1.trim.fq ,  $TMPDIR/${Sample}.R2.trim.fq using bwa-mem.."
    # bwa mem -t ${NSLOTS} -M $TMPDIR/BWAIndex/genome.fa $TMPDIR/${Sample}.R1.trim.fastq $TMPDIR/${Sample}.R2.trim.fastq | samtools view -bS - >  $TMPDIR/${Sample}.bam
#fi



#samtools rmdup
echo "Removing duplicates..."
out2=$(echo ${out1} | sed 's/\.bam$/.nodup.bam/')
echo ${out2}
picard MarkDuplicates INPUT=${out1} OUTPUT=${out2} METRICS_FILE="${out2}.dups.log" REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
#picard MarkDuplicates INPUT=Sample_N1.sorted.bam OUTPUT=Sample_N1.sorted.nodup.bam METRICS_FILE="sample.dups.log" REMOVE_DUPLICATES=true
# index
samtools index $out2


out2m=$(echo $out1 | sed 's/\.bam$/.nodup.noM.temp.bam/')

samtools idxstats $out2 | cut -f 1 | grep -v chrM | xargs samtools view -b $out2 > $out2m


#something odd happening
out2mb=$(echo $out1 | sed 's/\.bam$/.no.black.bam/')
#bedtools subtract -A -a $out2m -b $BLACK > $out2mb
# Remove multimapping and improper reads

out3=$(echo $out1 | sed 's/\.bam$/.nodup.noM.bam/')
samtools view -F 1804 -b ${out2m} > ${out3}
samtools index $out3



# histogram file
for w in 1000 500 200
do
    picard CollectInsertSizeMetrics I=$out3 O="${out3}.window${w}.hist_data" H="${out3}.window${w}.hist_graph.pdf" W=${w}
    done



## make bam with marked dups and generate PBC file for QC

samtools sort -n $p1 -o ${out1prefix}.nsort.bam

#sambamba sort --memory-limit 30GB \
#         --nthreads ${NSLOTS} --tmpdir ${TMPDIR} --out ${TMPDIR}/${Sample}/${Sample}.bam ${TMPDIR}/${Sample}/${Sample}.bam

samtools fixmate -r ${out1prefix}.nsort.bam ${out1prefix}.nsort.fixmate.bam
samtools view -F 1804 -f 2 -u  ${out1prefix}.nsort.fixmate.bam | samtools sort - > ${out1prefix}.filt.srt.bam


picard MarkDuplicates INPUT=${out1prefix}.filt.srt.bam OUTPUT=${out1prefix}.dupmark.bam METRICS_FILE=${out1prefix}.dup.qc VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false

samtools sort -n ${out1prefix}.dupmark.bam -o ${out1prefix}.srt.tmp.bam

dupmark_bam="${out1prefix}.srt.tmp.bam" 

PBC_QC="${out1prefix}.pbc.qc"

bedtools bamtobed -bedpe -i $dupmark_bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${PBC_QC}
