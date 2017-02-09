#!/bin/bash -l
#$ -N QCathena
#$ -j y
#$ -m e
#$ -cwd
### -l zenodotus=true
#$ -l athena=true
#$ -M ashley.doane@gmail.com
###$ -l h_rt=42:00:00
#$ -pe smp 2
#$ -l h_vmem=10G
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




cd ${path}
cd ${Sample}

ATHENA=1

if [ $ATHENA == 1 ]
then
    REFDIR="/athena/elementolab/scratch/asd2007/Reference"
    ANNOTDIR="/athena/elementolab/scratch/asd2007/Reference"
    #cp /home/asd2007/Scripts/picardmetrics.athena.conf ./picardmetrics.conf
else
    REFDIR="/zenodotus/dat01/melnick_bcell_scratch/asd2007/Reference"
    ANNOTDIR="/zenodotus/dat01/melnick_bcell_scratch/asd2007/Reference"
fi

#if [ $gtf_path == 1 ]
#then
    #gtf="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/mm10_UCSC_ref.gtf"
	REF="${REFDIR}/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex"
	REFbt2="${REFDIR}/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index"
  # REFbt="/home/asd20i07/dat02/asd2007/Reference/Homo_sapiens/UCSC/mm10/Sequence/BowtieIndex/genome"
  BLACK="${REFDIR}/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz" \
  #  BLACK="${REFDIR}/encodeBlack.bed"
    #BLACK="/athena/elementolab/scratch/asd2007/Reference/hg19/Anshul_Hg19UltraHighSignalArtifactRegions.bed"
    #chrsz="/athena/elementolab/scratch/asd2007/Reference/hg19.genome.chrom.sizes"
    #fetchChromSizes hg19 > hg19.chrom.sizes
    chrsz= $PWD/hg19.chrom.sizes
    RG="hg19"
    SPEC="hs"
    REFGen="/athena/elementolab/scratch/asd2007/bin/bcbio/genomes/Hsapiens/hg19/seq/"

#picardmetrics run -f "/home/asd2007/picardmetrics.athena.conf"  -o ${Sample}/QCmetrics/filtered ${Sample}/${Sample}.sorted.nodup.noM.bam
#cp ${Sample}/QCmetrics/filtered/*.EstimateLibraryComplexity.log ${Sample}/QCmetrics/${Sample}.picardcomplexity.qc
#cp ${Sample}/QCmetrics/filtered/*.EstimateLibraryComplexity.log ${Sample}/${Sample}.picardcomplexity.qc




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

python /home/asd2007/ATACseq/runFrip.py --workdir $PWD/${Sample} \
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














#Add file checker if fastq or gz

#copy the Samples with .gz extension to the tmp folder
#rsync -r -v -a -z $path/$file/*.fastq ./

#rsync -r -v -a -z $path/$file/*.fastq ./



#rsync -r -v -a -z $path/$file/* ./




#fastq-dump --split-files ${SRA}



#rsync -avP $TMPDIR/${Sample} $path/${Sample}
#rsync -r -a -v $TMPDIR/${Sample} $path/${Sample}


#############################




























# Don't forget to run `export -f module` first
#module add bedtools java picard-tools preseq python_anaconda r samtools ucsc_tools


# ~/dat02/asd2007/TOOLS/ataqc/run_ataqc.sh  /pbtech_mounts/oelab_store008/asd2007/dat02/asd2007/Projects/DataSets/atacData/atacEC3986/Sample_GCB_138/Sample_GCB_138 QCmetrics/atacqc Sample_GCB_138 Sample_GCB_138

# Directories and prefixes

#WORKDIR="/home/asd2007/dat02/asd2007/Projects/DataSets/atacData/atacEC3986/samples/Sample_GCB_138_FAST"
#OUTDIR="/home/asd2007/dat02/asd2007/Projects/DataSets/atacData/atacEC3986/samples/Sample_GCB_138_FAST"
##"/home/asd2007/dat02/asd2007/Projects/DataSets/atacData/atacEC3986/Sample_NB_138/Sample_NB_138/QCmetrics/atacqc"
#OUTPREFIX="Sampple_GCB_138_FAST"
##"QCmetrics"
#INPREFIX="Sample_GCB_FAST"
##"Sample_NB_138"
#GENOME='hg19' # This is the only genome that currently works
#
#
##SAMPLE=$(basename "$file")
##SAMPLE=${WORKDIR%%.*}
#
#echo "sample is" $SAMPLE
#
#SAMPLE=$INPREFIX
##"Sample_NB_138"
## Annotation files
#ANNOTDIR="/home/asd2007/melnick_bcell_scratch/asd2007/Reference"
#X="/home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/"
#DNASE_BED="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz"
#BLACKLIST_BED="${ANNOTDIR}/${GENOME}/Anshul_Hg19UltraHighSignalArtifactRegions.bed.gz"
#TSS_BED="${ANNOTDIR}/${GENOME}/hg19_RefSeq_stranded.bed.gz"
#REF_FASTA="${ANNOTDIR}/${GENOME}/encodeHg19Male.fa"
#PROM="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_prom_p2.bed.gz"
#ENH="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_enh_p2.bed.gz"
#REG2MAP="${ANNOTDIR}/${GENOME}/dnase_avgs_reg2map_p10_merged_named.pvals.gz"
#ROADMAP_META="${ANNOTDIR}/${GENOME}/eid_to_mnemonic.txt"
#
##WORKDIR="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/qc/rep2"
##OUTDIR="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/qc/rep2"
#outprefix="Sample_GCB_138_FAST"
#
##FASTQ1=$(ls ${WORKDIR}/../*R1.trim.fastq.gz)
#
#
##FASTQ2=$(ls ${WORKDIR}/../*R2.trim.fastq.gz)
#
#python $HOME/dat02/asd2007/TOOLS/ataqc/run_ataqc.py \
#        --workdir $WORKDIR \
#        --outdir $OUTDIR \
#        --outprefix $SAMPLE \
#        --genome hg19 \
#        --ref /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/encodeHg19Male/encodeHg19Male.fa \
#		    --tss /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/hg19_RefSeq_stranded.bed.gz \
#		    --dnase /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz \
#		    --blacklist /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz \
#		    --prom /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/reg2map_honeybadger2_dnase_prom_p2.bed.gz \
#		    --enh /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/reg2map_honeybadger2_dnase_enh_p2.bed.gz \
#		    --reg2map /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/dnase_avgs_reg2map_p10_merged_named.pvals.gz \
#		    --meta /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/eid_to_mnemonic.txt \
#		    --pbc "${WORKDIR}"/$SAMPLE*.pbc.qc \
#		    --fastq1 /home/asd2007/dat02/asd2007/Projects/DataSets/atacData/atacEC3986/samples/Sample_GCB_138_FAST/Sample_GCB_138_FAST.R1.trim.fastq.gz \
#        --fastq2 /home/asd2007/dat02/asd2007/Projects/DataSets/atacData/atacEC3986/samples/Sample_GCB_138_FAST/Sample_GCB_138_FAST.R2.trim.fastq.gz \
#		    --alignedbam /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/align/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.bam \
#		    --alignmentlog /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/qc/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.align.log \
#		    --coordsortbam /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/align/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.bam \
#		    --duplog /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/qc/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.dup.qc \
#		    --finalbam /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/align/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.nodup.bam \
#		    --finalbed /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/align/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.nodup.tn5.tagAlign.gz \
#		    --bigwig /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/signal/macs2/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.nodup.tn5.pf.pval.signal.bigwig \
#		    --peaks /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/macs2/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.nodup.tn5.pf.narrowPeak.gz \
#		    --naive_overlap_peaks /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/macs2/overlap/Sample_GCB_138.R1.trim.PE2SE.nodup.tn5_pooled.pf.pval0.1.500K.naive_overlap.narrowPeak.gz \
#		    --idr_peaks /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/idr/optimal_set/ppr.IDR0.1.filt.narrowPeak.gz

# SYS command. line 1112






: << COMMENTBLOCK







echo python /home/asd2007/ATACseq/run_ataqc.py --workdir $PWD/${Sample} \
    --outdir qc \
    --outprefix ${Sample} \
    --genome hg19 \
    --ref /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/encodeHg19Male/encodeHg19Male.fa \
    --tss /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/hg19_RefSeq_stranded.bed.gz \
    --dnase /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz \
    --blacklist /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz \
    --prom /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/reg2map_honeybadger2_dnase_prom_p2.bed.gz \
    --enh /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/reg2map_honeybadger2_dnase_enh_p2.bed.gz \
    --reg2map /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/dnase_avgs_reg2map_p10_merged_named.pvals.gz \
    --meta /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/eid_to_mnemonic.txt \
    --fastq1 ${F1} \
    --fastq2 ${F2} \
    --alignedbam ${ALIGNED_BAM} \
    --alignmentlog "${Sample}/${Sample}.align.log" \
    --coordsortbam "${Sample}/${Sample}.sorted.bam" \
    --duplog "${Sample}/${Sample}.dup.qc" \
    --pbc "${Sample}/${Sample}.pbc.qc" \
    --finalbam "${FINAL_BAM}" \
    --finalbed "${FINAL_BED}" \
    --bigwig "$Sample/$Sample.readsInPeaks.bin5.centered.smooth.150.max200.bw" \
    --peaks "${Sample}/${Sample}.tn5.pf.narrowPeak.gz" \
    --naive_overlap_peaks "${Sample}/pseudo_reps/${Sample}.nodup.tn5.pooled.pf.pval0.1.500K.naive_overlap.narrowPeak.gz" \
    --idr_peaks "${Sample}/pseudo_reps/${Sample}.nodup.tn5.pooled.pf.pval0.1.500K.naive_overlap.narrowPeak.gz"













                    --outprefix $OUTPREFIX \
                    --genome $GENOME \
                    --ref $REF_FASTA \
                    --alignmentlog $SAMPLE.bt2.log \
                    --tss $TSS_BED \
                    --dnase $DNASE_BED \
                    --blacklist $BLACKLIST_BED \
                    --prom $PROM \
                    --enh $ENH \
                    --reg2map $REG2MAP \
                    --meta $ROADMAP_META \
                    --inprefix $INPREFIX \
                    --finalbam $SAMPLE.subsample.sorted.nodup.noM.black.bam \
                    --finalbed $SAMPLE.tag.narrow_peaks.narrowPeak \
                    --coordsortbam $SAMPLE.subsample.sorted.nodup.noM.black.bam \
                    --alignedbam $SAMPLE.subsample.sorted.bam \
                    --duplog $SAMPLE.sorted.duplicate_metrics \
                    --pbc $SAMPLE.subsample.pbc.qc \
                    --bigwig $SAMPLE.smooth151.center.extend.fpkm.bw \
                    --fastq1 $FASTQ1 \
                    --fastq2 $FASTQ2 \

                    --peaks $SAMPLE.narrow_peaks.narrowPeak 


 /home/asd2007/melnick_bcell_scratch/asd2007/bin/bds_atac/ataqc/run_ataqc.py \
		    --workdir /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/qc/rep2 \
		    --outdir /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/qc/rep2 \
		    --outprefix Sample_GCB_138_FAST.R1.trim.PE2SE \
		    --genome hg19 \
		    --ref /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/encodeHg19Male/encodeHg19Male.fa \
		    --tss /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/hg19_RefSeq_stranded.bed.gz \
		    --dnase /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz \
		    --blacklist /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz \
		    --prom /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/reg2map_honeybadger2_dnase_prom_p2.bed.gz \
		    --enh /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/reg2map_honeybadger2_dnase_enh_p2.bed.gz \
		    --reg2map /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/dnase_avgs_reg2map_p10_merged_named.pvals.gz \
		    --meta /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/eid_to_mnemonic.txt \
		    --pbc /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/qc/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.nodup.pbc.qc\
		    --fastq1 /home/asd2007/dat02/asd2007/Projects/DataSets/atacData/atacEC3986/samples/Sample_GCB_138_FAST/Sample_GCB_138_FAST.R1.trim.fastq.gz --fastq2 /home/asd2007/dat02/asd2007/Projects/DataSets/atacData/atacEC3986/samples/Sample_GCB_138_FAST/Sample_GCB_138_FAST.R2.trim.fastq.gz \
		    --alignedbam /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/align/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.bam \
		    --alignmentlog /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/qc/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.align.log \
		    --coordsortbam /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/align/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.bam \
		    --duplog /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/qc/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.dup.qc \
		    --finalbam /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/align/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.nodup.bam \
		    --finalbed /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/align/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.nodup.tn5.tagAlign.gz \
		    --bigwig /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/signal/macs2/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.nodup.tn5.pf.pval.signal.bigwig \
		    --peaks /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/macs2/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.nodup.tn5.pf.narrowPeak.gz \
		    --naive_overlap_peaks /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/macs2/overlap/Sample_GCB_138.R1.trim.PE2SE.nodup.tn5_pooled.pf.pval0.1.500K.naive_overlap.narrowPeak.gz \
		    --idr_peaks /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/idr/optimal_set/ppr.IDR0.1.filt.narrowPeak.gz \
		    

# SYS command. line 1112




    # Mode 2: Define every possible QC file
    parser.add_argument('--fastq1',
                        help='First set of reads if paired end, \
                              or the single end reads')
    parser.add_argument('--fastq2',
                        help='Second set of reads if paired end')
    parser.add_argument('--alignedbam', help='BAM file from the aligner')
    parser.add_argument('--alignmentlog', help='Alignment log')
    parser.add_argument('--coordsortbam', help='BAM file sorted by coordinate')
    parser.add_argument('--duplog', help='Picard duplicate metrics file')
    parser.add_argument('--pbc', help='ENCODE library complexity metrics file')
    parser.add_argument('--finalbam', help='Final filtered BAM file')
    parser.add_argument('--finalbed',
                        help='Final filtered alignments in BED format')
    parser.add_argument('--bigwig',
                        help='Final bigwig')
    parser.add_argument('--peaks',
                        help='Peak file')
    parser.add_argument('--naive_overlap_peaks',
                        default=None, help='Naive overlap peak file')
    parser.add_argument('--idr_peaks',
                        default=None, help='IDR peak file')
    parser.add_argument('--use_sambamba_markdup', action='store_true',
                        help='Use sambamba markdup instead of Picard')

COMMENTBLOCK
