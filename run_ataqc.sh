#!/bin/bash -l

# This script is an example test of the ataqc package

# Don't forget to run `export -f module` first
#module add bedtools java picard-tools preseq python_anaconda r samtools ucsc_tools


# ~/dat02/asd2007/TOOLS/ataqc/run_ataqc.sh  /pbtech_mounts/oelab_store008/asd2007/dat02/asd2007/Projects/DataSets/atacData/atacEC3986/Sample_GCB_138/Sample_GCB_138 QCmetrics/atacqc Sample_GCB_138 Sample_GCB_138

# Directories and prefixes
WORKDIR="/home/asd2007/dat02/asd2007/Projects/DataSets/atacData/atacEC3986/samples/Sample_GCB_138_FAST"
OUTDIR="/home/asd2007/dat02/asd2007/Projects/DataSets/atacData/atacEC3986/samples/Sample_GCB_138_FAST"
#"/home/asd2007/dat02/asd2007/Projects/DataSets/atacData/atacEC3986/Sample_NB_138/Sample_NB_138/QCmetrics/atacqc"
OUTPREFIX="Sampple_GCB_138_FAST"
#"QCmetrics"
INPREFIX="Sample_GCB_FAST"
#"Sample_NB_138"
GENOME='hg19' # This is the only genome that currently works


#SAMPLE=$(basename "$file")
#SAMPLE=${WORKDIR%%.*}

echo "sample is" $SAMPLE

SAMPLE=$INPREFIX
#"Sample_NB_138"
# Annotation files
ANNOTDIR="/home/asd2007/melnick_bcell_scratch/asd2007/Reference"
X="/home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/"
DNASE_BED="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz"
BLACKLIST_BED="${ANNOTDIR}/${GENOME}/Anshul_Hg19UltraHighSignalArtifactRegions.bed.gz"
TSS_BED="${ANNOTDIR}/${GENOME}/hg19_RefSeq_stranded.bed.gz"
REF_FASTA="${ANNOTDIR}/${GENOME}/encodeHg19Male.fa"
PROM="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_prom_p2.bed.gz"
ENH="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_enh_p2.bed.gz"
REG2MAP="${ANNOTDIR}/${GENOME}/dnase_avgs_reg2map_p10_merged_named.pvals.gz"
ROADMAP_META="${ANNOTDIR}/${GENOME}/eid_to_mnemonic.txt"

#WORKDIR="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/qc/rep2"
#OUTDIR="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/qc/rep2"
outprefix="Sample_GCB_138_FAST"

#FASTQ1=$(ls ${WORKDIR}/../*R1.trim.fastq.gz)


#FASTQ2=$(ls ${WORKDIR}/../*R2.trim.fastq.gz)

python $HOME/dat02/asd2007/TOOLS/ataqc/run_ataqc.py \
        --workdir $WORKDIR \
        --outdir $OUTDIR \
        --outprefix $SAMPLE \
        --genome hg19 \
        --ref /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/encodeHg19Male/encodeHg19Male.fa \
		    --tss /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/hg19_RefSeq_stranded.bed.gz \
		    --dnase /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz \
		    --blacklist /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz \
		    --prom /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/reg2map_honeybadger2_dnase_prom_p2.bed.gz \
		    --enh /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/reg2map_honeybadger2_dnase_enh_p2.bed.gz \
		    --reg2map /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/dnase_avgs_reg2map_p10_merged_named.pvals.gz \
		    --meta /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/eid_to_mnemonic.txt \
		    --pbc "${WORKDIR}"/$SAMPLE*.pbc.qc \
		    --fastq1 /home/asd2007/dat02/asd2007/Projects/DataSets/atacData/atacEC3986/samples/Sample_GCB_138_FAST/Sample_GCB_138_FAST.R1.trim.fastq.gz \
        --fastq2 /home/asd2007/dat02/asd2007/Projects/DataSets/atacData/atacEC3986/samples/Sample_GCB_138_FAST/Sample_GCB_138_FAST.R2.trim.fastq.gz \
		    --alignedbam /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/align/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.bam \
		    --alignmentlog /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/qc/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.align.log \
		    --coordsortbam /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/align/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.bam \
		    --duplog /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/qc/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.dup.qc \
		    --finalbam /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/align/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.nodup.bam \
		    --finalbed /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/align/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.nodup.tn5.tagAlign.gz \
		    --bigwig /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/signal/macs2/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.nodup.tn5.pf.pval.signal.bigwig \
		    --peaks /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/macs2/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.nodup.tn5.pf.narrowPeak.gz \
		    --naive_overlap_peaks /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/macs2/overlap/Sample_GCB_138.R1.trim.PE2SE.nodup.tn5_pooled.pf.pval0.1.500K.naive_overlap.narrowPeak.gz \
		    --idr_peaks /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/idr/optimal_set/ppr.IDR0.1.filt.narrowPeak.gz

# SYS command. line 1112






: << COMMENTBLOCK

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
