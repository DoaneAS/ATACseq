#!/bin/bash

#Sample_TH29_3.narrow.p0.1_peaks.narrowPeak


# ========================
# Create pseudoReplicates
# =======================
FINAL_BEDPE_FILE=$1
SAMPLE=$2
SPEC=$3

chrsz="/home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19.genome.chrom.sizes"
# Get total number of read pairs
nlines=$( zcat ${FINAL_BEDPE_FILE} | wc -l  )
nlines=$(( (nlines + 1) / 2  ))

# Shuffle and split BEDPE file into 2 equal parts

# ========================
# Create pseudoReplicates
# =======================



mkdir -p ${SAMPLE}/pseudo_reps


zcat $FINAL_BEDPE_FILE | shuf | split -d -l $((nlines)) - $SAMPLE/pseudo_reps/$SAMPLE.nodup.

# SYS command. line 303

awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' "$SAMPLE/pseudo_reps/$SAMPLE.nodup.00" | \
    gzip -c > $SAMPLE/pseudo_reps/$SAMPLE.nodup.pr1.tagAlign.gz

rm -f "$TMPDIR/$SAMPLE/pseudo_reps/$SAMPLE.nodup.00"


awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' "$SAMPLE/pseudo_reps/$SAMPLE.nodup.01" | \
		gzip -c > gzip -c > $SAMPLE/pseudo_reps/$SAMPLE.nodup.pr2.tagAlign.gz


rm -f "$TMPDIR/$SAMPLE/pseudo_reps/$SAMPLE.nodup.01"
# SYS command. line 308

POOLED_PREFIX="$SAMPLE/$SAMPLE.nodup"
POOLED_TA_FILE="$SAMPLE/$SAMPLE.nodup.tagAlign.gz"

PR_PREFIX="$SAMPLE/pseudo_reps/$SAMPLE.nodup"
PR1_TA_FILE="$SAMPLE/pseudo_reps/$SAMPLE.nodup.pr1.tagAlign.gz"
PR2_TA_FILE="$SAMPLE/pseudo_reps/$SAMPLE.nodup.pr2.tagAlign.gz"


zcat ${POOLED_TA_FILE} |  awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' | gzip -c > ${POOLED_PREFIX}.pooled.tn5.tagAlign.gz
zcat ${PR1_TA_FILE} |  awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' | gzip -c > ${PR_PREFIX}.pr1.tn5.tagAlign.gz
zcat ${PR2_TA_FILE} |  awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' | gzip -c > ${PR_PREFIX}.pr2.tn5.tagAlign.gz
# Get total number of read pairs

# Shuffle and split BEDPE file into 2 equal parts


macs2 callpeak -t   ${POOLED_PREFIX}.pooled.tn5.tagAlign.gz -f BED -n  ${POOLED_PREFIX}.pooled -g $SPEC  --nomodel --shift -75 --extsize 150 --keep-dup all --bdg --SPMR --call-summits -p 1e-1
sort -k 8gr,8gr ${POOLED_PREFIX}.pooled_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | gzip -nc > ${POOLED_PREFIX}.pooled.tn5.narrowPeak.gz

macs2 callpeak -t   ${PR_PREFIX}.pr1.tn5.tagAlign.gz -f BED -n  ${PR_PREFIX}.pr1 -g $SPEC  --nomodel --shift -75 --extsize 150 --keep-dup all --bdg --SPMR --call-summits -p 1e-1
sort -k 8gr,8gr ${PR_PREFIX}.pr1_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | gzip -nc > ${PR_PREFIX}.pr1.tn5.narrowPeak.gz

macs2 callpeak -t ${PR_PREFIX}.pr2.tn5.tagAlign.gz -f BED -n  ${PR_PREFIX}.pr2 -g $SPEC  --nomodel --shift -75 --extsize 150 --keep-dup all --bdg --SPMR --call-summits -p 1e-1


sort -k 8gr,8gr ${PR_PREFIX}.pr2_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | gzip -nc > ${PR_PREFIX}.pr2.tn5.narrowPeak.gz

#prefix=Sample_TH29_3
pval_thresh=0.01




##out
#Narrowpeak file ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.narrowPeak.gz
#Broadpeak file ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.broadPeak.gz
#Gappedpeak file ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.gappedPeak.gz
#Fold-enrichment bigWig file ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.fc.signal.bw
#-log10(pvalue) bigWig file ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.pval.signal.bw
#


zcat ${POOLED_PREFIX}.pooled.tn5.narrowPeak.gz | sort -grk8 | head -n 500000 | gzip -c > ${POOLED_PREFIX}.pooled.tn5.pval0.1.500k.narrowPeak.gz

zcat ${PR_PREFIX}.pr1.tn5.narrowPeak.gz | sort -grk8 | head -n 500000 | gzip -c > ${PR_PREFIX}.pr1.tn5.pval0.1.500k.narrowPeak.gz


zcat ${PR_PREFIX}.pr2.tn5.narrowPeak.gz | sort -grk8 | head -n 500000 | gzip -c > ${PR_PREFIX}.pr2.tn5.pval0.1.500k.narrowPeak.gz

POOLED=${POOLED_PREFIX}.pooled.tn5.pval0.1.500k.narrowPeak.gz


#zcat ${SAMPLE}/${SAMPLE}.tn5.pf.narrowPeak.gz  | sort -grk8 | head -n 500000 | gzip -c > ${PR_PREFIX}.pooled.tn5.pval0.1.500k.narrowPeak.gz



intersectBed -wo -a <(zcat -f ${POOLED}) -b <(zcat -f ${PR_PREFIX}.pr1.tn5.pval0.1.500k.narrowPeak.gz) | awk 'BEGIN{FS="    ";OFS="    "} {s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq | intersectBed -wo -a stdin -b <(zcat -f ${PR_PREFIX}.pr2.tn5.pval0.1.500k.narrowPeak.gz) | awk 'BEGIN{FS="    ";OFS="    "} {s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq | gzip -c > ${PR_PREFIX}.tn5.pooled.pf.pval0.1.500K.naive_overlap.narrowPeak.gz

# SYS command. line 109


# SYS command. line 112


# SYS command. line 114


# SYS command. line 116

#TASKTIME=$[$(date +%s)-${STARTTIME}]; echo "Task has finished (${TASKTIME} seconds)."
