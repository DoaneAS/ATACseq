#!/bin/bash
# read from command line which files to align
p1=$1

BLACK="/home/asd2007/melnick_bcell_scratch/asd2007/Reference/encodeBlack.bed"
# help
if [ -z "$p1"  ]
then
    echo "This will sort bam, remove mitochondria/duplicates and histogram insert size"
    echo "$(basename $0) <bamFile>"
    echo "<samFile>: Required bam input file"
    exit
fi

# convert to bam and sort
echo "Sorting..."
out1prefix=$(echo $p1 | sed 's/\.bam$//')
out1="${out1prefix}.sorted.bam"
echo ${out1}
samtools view -S -u $p1 | samtools sort - > ${out1}
#index
samtools index $out1

#samtools rmdup
echo "Removing duplicates..."
out2=$(echo ${out1} | sed 's/\.bam$/.nodup.bam/')
echo ${out2}
picard MarkDuplicates INPUT=${out1} OUTPUT=${out2} METRICS_FILE="${out2}.dups.log" REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
#picard MarkDuplicates INPUT=Sample_N1.sorted.bam OUTPUT=Sample_N1.sorted.nodup.bam METRICS_FILE="sample.dups.log" REMOVE_DUPLICATES=true
# index
samtools index $out2


out2m=$(echo $out1 | sed 's/\.bam$/.nodup.noM.bam/')

samtools idxstats $out2 | cut -f 1 | grep -v chrM | xargs samtools view -b $out2 > $out2m


#something odd happening
out2mb=$(echo $out1 | sed 's/\.bam$/.no.black.bam/')
bedtools subtract -A -a $out2m -b $BLACK > $out2mb
# Remove multimapping and improper reads

out3=$(echo $out1 | sed 's/\.bam$/.nodup.noM.black.bam/')
samtools view -F 1804 -b ${out2mb} > ${out3}
samtools index $out3



# histogram file
for w in 1000 500 200
do
    picard CollectInsertSizeMetrics I=$out3 O="${out3}.window${w}.hist_data" H="${out3}.window${w}.hist_graph.pdf" W=${w}
    done
