

#source activate bds_atac_py3



REP1_PEAK_FILE=$1
REP2_PEAK_FILE=$2

POOLED_PEAK_FILE=$3

IDR_OUTPUT=$4
IDR_THRESH=0.1

BLACKLIST="/athena/elementolab/scratch/asd2007/Reference/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz"

#zcat ${REP1_TA_FILE} ${REP2_TA_FILE} | gzip -c > ${POOLED_TA_FILE}


# =============================
# Perform IDR analysis.
# Generate a plot and IDR output with additional columns including IDR scores.
# =============================
idr --samples ${REP1_PEAK_FILE} ${REP2_PEAK_FILE} --peak-list ${POOLED_PEAK_FILE} --input-file-type narrowPeak --output-file ${IDR_OUTPUT} --rank p.value --soft-idr-threshold ${IDR_THRESH} --plot --use-best-multisummit-IDR

# =============================
# Get peaks passing IDR threshold of 10%
# =============================
IDR_THRESH_TRANSFORMED=$(awk -v p=${IDR_THRESH} 'BEGIN{print -log(p)/log(10)}')

awk 'BEGIN{OFS="\t"} $12>='"${IDR_THRESH_TRANSFORMED}"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' ${IDR_OUTPUT} | sort | uniq | sort -k7n,7n | gzip -c > ${IDR_OUTPUT}.IDR0.1.narrowPeak.gz

NPEAKS_IDR=$(zcat ${IDR_OUTPUT}.IDR0.1.narrowPeak.gz | wc -l)

cat $NPEAKS_IDR > ${IDR_OUTPUT}.numIDR.txt

# =============================
# Filter using black list
# =============================
bedtools intersect -v -a ${IDR_OUTPUT}.IDR0.1.narrowPeak.gz -b ${BLACKLIST} | gzip -c > ${IDR_OUTPUT}.IDR0.1.filt.narrowPeak.gz



source deactivate

