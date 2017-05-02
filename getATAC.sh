#!/bin/bash -l
#$ -N ATACathena
#$ -j y
#$ -m e
#$ -cwd
###$ -l zeno=true
#$ -l athena=true
#$ -M ashley.doane@gmail.com
#$ -l h_rt=48:00:00
#$ -pe smp 4-8
#$ -l h_vmem=6G
#$ -R y
#$ -o /home/asd2007/joblogs

########### SETTINGS ###########
TRIM="YES" # 
NUC=0 # run nucleoatac, extended runtime required
BT2ALN="YES" # bt2 recommnded!
ATHENA=1
##############################

while [[ $# -gt 1 ]]
do
    key="$1"

    case $key in
        -f|--folderpath)
            FOLDERPATH="$2"
            shift # past argument
            ;;
        -g|--genome)
            GENOME="$2"
            shift # past argument
            ;;
        -t|--trim)
            TRIM="$2"
            shift # past argument
            ;;
        --align)
            BT2ALN=YES
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
echo BT2 ALN    = "${BT2ALN}"
if [[ -n $1 ]]; then
    echo "Last line of FILEPATH specified as non-opt/last argument:"
    tail -1 $1
fi


BT2LAN=YES

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

#pigz -p $NSLOTS -c *.R1.trim.fastq > $TMPDIR/${Sample}.R1.trim.fastq.gz
#pigz -p $NSLOTS -c *.R2.trim.fastq > $TMPDIR/${Sample}.R2.trim.fastq.gz
#rm *trim.fastq

#parallel -j ${NSLOTS} 'fastq-dump --split-FILEPATHs {}' ::: *.sra

echo "Processing  $Sample ..."

#Figuring out the reference genome


if [[ $ATHENA == 1 ]] ; then
    REFDIR="/athena/elementolab/scratch/asd2007/reference"
    ANNOTDIR="/athena/elementolab/scratch/asd2007/reference"
    #PICARDCONF="/home/asd2007/Scripts/picardmetrics.athena.conf"
    export PATH="/home/asd2007/anaconda2/bin:$PATH"
else
    REFDIR="/zenodotus/dat01/melnick_bcell_scratch/asd2007/reference"
    ANNOTDIR="/zenodotus/dat01/melnick_bcell_scratch/asd2007/reference"
    PICARDCONF="/home/asd2007/Scripts/picardmetrics.conf"
fi

# get genome args
if [[ $GENOME == "hg19" ]] ; then
	  REF="${REFDIR}/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex"
	  REFbt2="${REFDIR}/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index"
    BLACK="${REFDIR}/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
    #fetchChromSizes hg19 > hg19.chrom.sizes
    chrsz="/athena/elementolab/scratch/asd2007/reference/hg19.chrom.sizes"
    cp $chrsz $PWD/hg19.chrom.sizes
    RG="hg19"
    SPEC="hs"
    REFGen="/athena/elementolab/scratch/asd2007/local/share/bcbio/genomes/Hsapiens/hg19/seq/hg19.fa"
elif [[ $GENOME == "hg38" ]] ; then
    echo "genome is ${GENOME}"
    DNASE_BED="${ANNOTDIR}/${GENOME}/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz"
    BLACK="/athena/elementolab/scratch/asd2007/reference/hg38/hg38.blacklist.bed.gz"
    PICARDCONF="/athena/elementolab/scratch/asd2007/reference/hg38/picardmetrics.conf"
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


#rsync -avP ${REF} $TMPDIR/#
#rsync -av ${REFbt2} $TMPDIR/

mkdir -p ${TMPDIR}/${Sample}

echo "ls of pwd is"
ls -lrth


#find *_L00*_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz$//' |parallel 'trimAdapters.py -a {}_R1_001.fastq.gz -b {}_R2_001.fastq.gz'

#'cutadapt -a adaptors_to_trim -A adaptors_to_trim -q 20 --minimum-length 5 -o {}_R1_cutadapt.fastq.gz -p {}_R2_cutadapt.fastq.gz {}_R1.fastq.gz {}_R2.fastq.gz &> {}.cutadapt'

#cat *R1.trim.fastq > ${Sample}.R1.trim.fastq
#cat *R2.trim.fastq > ${Sample}.R2.trim.fastq

if [ -f "${TMPDIR}/${Sample}.R1.trim.fastq.gz" ];
then
    TRIM=NO
else
    TRIM=YES
fi



if [[ $TRIM == YES ]]; then
    echo "Trimming adapter sequences, with command..."
    #find *_L00*_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz$//' |parallel -j ${NSLOTS} 'trimAdapters.py -a {}_R1_001.fastq.gz -b {}_R2_001.fastq.gz'
    find *_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz$//' |parallel -j ${NSLOTS} 'trimAdapters.py -a {}_R1_001.fastq.gz -b {}_R2_001.fastq.gz'
    cat *_R1_001.trim.fastq > ${Sample}.R1.trim.fastq
    cat *_R2_001.trim.fastq > ${Sample}.R2.trim.fastq
    echo "completed trimming"
    pigz -p $NSLOTS -c $TMPDIR/${Sample}.R1.trim.fastq > $TMPDIR/${Sample}.R1.trim.fastq.gz
    pigz -p $NSLOTS -c $TMPDIR/${Sample}.R2.trim.fastq > $TMPDIR/${Sample}.R2.trim.fastq.gz
    rsync -av $TMPDIR/${Sample}.R1.trim.fastq.gz ${FOLDERPATH}/${Sample}
    rsync -av $TMPDIR/${Sample}.R2.trim.fastq.gz ${FOLDERPATH}/${Sample}
else
    echo "will not perform adapter sequence trimming of reads"
    #gunzip *.gz
      #cat *R1.trim.fastq.gz >  $TMPDIR/${Sample}.R1.trim.fastq.gz
      #cat *R2.trim.fastq.gz > $TMPDIR/${Sample}.R2.trim.fastq.gz
      #mkdir R1
      #mkdir R2
      #cat ${Sample}.R1.trim.fastq.gz >  $TMPDIR/${Sample}.R1.trim.fastq.gz
      #cat ${Sample}.R2.trim.fastq.gz > $TMPDIR/${Sample}.R2.trim.fastq.gz
      #cho "fastq FILEPATHs to align"
      #ls $TMPDIR/${Sample}.R*.trim.fastq.gz
fi




if [[ $BT2ALN == YES ]] ; then
    echo "----------bowtie2 aligning-------------"
    bowtie2 -X 2000 --threads ${NSLOTS} -x ${bwt2_idx}  \
           -1 $TMPDIR/${Sample}.R1.trim.fastq.gz -2 $TMPDIR/${Sample}.R2.trim.fastq.gz 2> ${Sample}/${Sample}.align.log | samtools view -bS - > $TMPDIR/${Sample}/${Sample}.bam
    echo $(cat ${Sample}/${Sample}.align.log)
else
    echo "----------bwa-mem aligning-------------"
    # echo "aligning : $TMPDIR/${Sample}.R1.trim.fq ,  $TMPDIR/${Sample}.R2.trim.fq using bwa-mem.."
   # bwa mem -t ${NSLOTS} -M ${REF} $TMPDIR/${Sample}.R1.trim.fastq.gz $TMPDIR/${Sample}.R2.trim.fastq.gz | samtools view -bS - >  $TMPDIR/${Sample}/${Sample}.bam
fi

##bwa mem -t ${NSLOTS} -M ${REF} $TMPDIR/${Sample}.R1.trim.fastq.gz $TMPDIR/${Sample}.R2.trim.fastq.gz | samtools view -bS - >  $TMPDIR/${Sample}/${Sample}.bam
echo "----------Processing alignment and filtering for duplicates and mitochondrial mapping reads-----------------"


spack load bedtools


rsync -a -v $TMPDIR/${Sample} $FOLDERPATH/${Sample}

processBamAlignment.sh $TMPDIR/${Sample}/${Sample}.bam ${BLACK}

sambamba sort --memory-limit 30GB -n \
         -t ${NSLOTS} --out $TMPDIR/${Sample}/${Sample}.nsorted.nodup.noM.bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam

# ALL bams below have blackList regions removed

FINAL_BAM='${TMPDIR}/${Sample}/${Sample}.sorted.nodup.noM.black.bam'


samtools fixmate $TMPDIR/${Sample}/${Sample}.nsorted.nodup.noM.bam $TMPDIR/${Sample}/${Sample}.nsorted.fixmate.nodup.noM.bam

convertBAMtoBED.sh $TMPDIR/${Sample}/${Sample}.nsorted.fixmate.nodup.noM.bam

samtools view -H $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.bam | grep chr | grep -v chrM | /home/ole2001/PERL_SCRIPTS/columns.pl 1 2 | sed 's/SN://' | sed 's/LN://' > chrom.sizes


echo "----------Finsihed Processing Bam Files-----------------"

#rsync -avP $TMPDIR/${Sample} $FOLDERPATH/${Sample}
rsync -r -a -v $TMPDIR/${Sample} $FOLDERPATH/${Sample}


echo "------------------------------------------Call Peaks with MACS2--------------------------------------------"

cp $TMPDIR/${Sample}/${Sample}.nsorted.fixmate.nodup.noM.tn5.tagAlign.gz  $TMPDIR/${Sample}/${Sample}.nodup.tn5.tagAlign.gz

cp $TMPDIR/${Sample}/${Sample}.nsorted.fixmate.nodup.noM.bedpe.gz $TMPDIR/${Sample}/${Sample}.nodup.bedpe.gz

macs2 callpeak -t  $TMPDIR/${Sample}/${Sample}.nodup.tn5.tagAlign.gz -f BED -n $TMPDIR/${Sample}/${Sample}.tag.narrow -g $SPEC  --nomodel --shift -75 --extsize 150 --keep-dup all --bdg --SPMR --call-summits -p 1e-2

/home/asd2007/ATACseq/narrowpeak.py $TMPDIR/${Sample}/${Sample}.tag.narrow_peaks.narrowPeak $TMPDIR/${Sample}/${Sample}.tn5.narrowPeak.gz 

sort -k 8gr,8gr $TMPDIR/${Sample}/${Sample}.tag.narrow_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | gzip -nc > $TMPDIR/${Sample}/${Sample}.tn5.pf.narrowPeak.gz

zcat $TMPDIR/${Sample}/${Sample}.tn5.pf.narrowPeak.gz | sort -grk8 | head -n 500000 | gzip -c > $$TMPDIR/${Sample}/${Sample}.tn5.pval0.01.500k.narrowPeak.gz

## BROAD peaks ##
macs2 callpeak -t  $TMPDIR/${Sample}/${Sample}.nsorted.fixmate.nodup.noM.tn5.tagAlign.gz -f BED -n $TMPDIR/${Sample}/${Sample}.tag.broad -g $SPEC  --nomodel --shift -75 --extsize 150 --keep-dup all --broad --broad-cutoff 0.1

sort -k 8gr,8gr $TMPDIR/${Sample}/${Sample}.tag.broad_peaks.broadPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | gzip -nc > $TMPDIR/${Sample}/${Sample}.tn5.pf.broadPeak.gz

sort -k 14gr,14gr $TMPDIR/${Sample}/${Sample}.tag.broad_peaks.gappedPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | gzip -nc > $TMPDIR/${Sample}/${Sample}.tn5.pf.gappedPeak.gz

mkdir -p $TMPDIR/${Sample}/peaks


cp $TMPDIR/${Sample}/${Sample}.tag*Peak  $TMPDIR/${Sample}/peaks/
cp $TMPDIR/${Sample}/${Sample}.tn5*Peak.gz  $TMPDIR/${Sample}/peaks/




getFrip.sh ${TMPDIR}/${Sample}/${Sample}.sorted.nodup.noM.black.bam $TMPDIR/${Sample}/${Sample}.tag.broad_peaks.broadPeak




rsync -r -a -v $TMPDIR/${Sample} $FOLDERPATH/${Sample}
echo "------------------------------------------Call Peaks with MACS2 ON PSEUDOREPLICTES--------------------------------------------"

callPeaks.sh $TMPDIR/${Sample}/${Sample}.nodup.bedpe.gz ${Sample} $SPEC


rsync -r -a -v $TMPDIR/${Sample} $FOLDERPATH/${Sample}


mkdir -p $TMPDIR/${Sample}/bamPE


macs2 callpeak -t ${TMPDIR}/${Sample}/${Sample}.sorted.nodup.noM.black.bam -f BAMPE -n $TMPDIR/${Sample}/bamPE/${Sample}.bamPE -g ${SPEC} --nomodel --keep-dup all --broad --broad-cutoff 0.1  --bdg --SPMR 


## get FE over background
~/ATACseq/getFC.sh $TMPDIR/${Sample}/bamPE/${Sample}.bamPE ${chrsz}


rsync -r -a -v $TMPDIR/${Sample}/bamPE $FOLDERPATH/${Sample}/


prefix=$TMPDIR/${Sample}/${Sample}.narrow


echo "----------- compute QC stats  ------------------"

mkdir -p ${Sample}/QCmetrics

mkdir -p ${Sample}/QCmetrics/raw

mkdir -p ${Sample}/QCmetrics/filtered

pyatac sizes --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.bam --upper 1000 --out $TMPDIR/${Sample}/QCmetrics/${Sample}.filtered


picardmetrics run -f $PICARDCONF -o $TMPDIR/${Sample}/QCmetrics/filtered $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam

picardmetrics run -f $PICARDCONF -o $TMPDIR/${Sample}/QCmetrics/raw $TMPDIR/${Sample}/${Sample}.sorted.bam

cp $TMPDIR/${Sample}/QCmetrics/raw/*.EstimateLibraryComplexity.log $TMPDIR/${Sample}/QCmetrics/${Sample}.picardcomplexity.qc
#cp $TMPDIR/${Sample}/QCmetrics/filtered/*.EstimateLibraryComplexity.log $TMPDIR/${Sample}/${Sample}.picardcomplexity.qc

samtools flagstat  $TMPDIR/${Sample}/${Sample}.sorted.bam > $TMPDIR/${Sample}/QCmetrics/raw/${Sample}.sorted.flagstat.txt
samtools flagstat $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.bam > $TMPDIR/${Sample}/QCmetrics/filtered/${Sample}.sorted.nodup.noM.flagstat.txt

getchrM $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.bam > $TMPDIR/${Sample/}QCmetrics/${Sample}.sorted.nodup.noM.bam.chrM.txt
getchrM $TMPDIR/${Sample}/${Sample}.sorted.bam > $TMPDIR/${Sample}/QCmetrics/${Sample}.sorted.bam.chrM.txt


echo "----------- compute fragment midpoint coverage per bp ------------------"
#pyatac cov --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.bam --out $TMPDIR/${Sample}/${Sample}.pyatac.cov --cores ${NSLOTS}

#pyatac cov --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.bam --out $TMPDIR/${Sample}/${Sample}.pyatac.125bp.cov --cores ${NSLOTS} --upper 125

#bamCoverage --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.bam --binSize 5 \
#    --outFileFormat bigwig --smoothLength 150 \
#    --normalizeUsingRPKM \
#    --maxFragmentLength 300 \
#    -o $TMPDIR/${Sample}/${Sample}.smooth151.center.extend.fpkm.max300.bw --centerReads --extendReads --numberOfProcessors ${NSLOTS}


bamCoverage --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.bam --binSize 5 \
            --outFileFormat bigwig --smoothLength 150 \
            --normalizeUsingRPKM \
            --maxFragmentLength 150 \
            -o $TMPDIR/${Sample}/${Sample}.smooth150.center.extend.fpkm.max150.bw --centerReads --extendReads --numberOfProcessors ${NSLOTS}


#bamCoverage --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.bam --binSize 20 \
#    --outFileFormat bigwig --smoothLength 150 \
#    --normalizeUsingRPKM \
#    -o $TMPDIR/${Sample}/${Sample}.bin20.smooth150.fpkm.bw --numberOfProcessors ${NSLOTS}

rsync -r -a -v $TMPDIR/${Sample} $FOLDERPATH/${Sample}




# Peak atlas using Chang data and B cells
#BED="/athena/elementolab/scratch/asd2007/Projects/DataSets/atacData/ATAC-AML/AML/cellAtlas.bed"

#BED="/athena/elementolab/scratch/asd2007/Projects/HOMER/Homer.input/BCellAtlas.atacSignal_sorted.bed"


#pyatac counts --bam  $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.bam --bed ${BED} --out $TMPDIR/${Sample}/${Sample}.150.ins --upper 150
##bedtools multicov -split -bams $BAMS -bed $BED

#pyatac counts --bam  $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam --bed ${BED} --out $TMPDIR/${Sample}/${Sample}.ins

#LIBS=$(zcat ${Sample}/${Sample}.ins.counts.txt.gz  | awk '{sum +=$1} END  {printf sum}')

#let KM=2000000 #1 million * bin size in kb
#let LIBZ=$LIBS/$KM

#bamCoverage --bam  $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.bam --binSize 20 \
#    --outFileFormat bigwig --smoothLength 150  \
#    --scaleFactor $LIBZ \
#    -o ${TMPDIR}/${Sample}/$Sample.readsInPeaks.bin20.centered.smooth.150.max150f.bw --numberOfProcessors ${NSLOTS} \
#    --maxFragmentLength 150 \
#    --centerReads --extendReads



#bamCoverage --bam ${TMPDIR}/${Sample}/${Sample}.sorted.nodup.noM.bam --binSize 20 \
#    --outFileFormat bigwig --smoothLength 150  \
#    --scaleFactor $LIBZ \
#    -o ${TMPDIR}/$Sample/$Sample.readsInPeaks.bin20.centered.smooth.150.max500f.bw --numberOfProcessors ${NSLOTS} \
#    --maxFragmentLength 500 \
#    --centerReads --extendReads


let KM=50000 #1 million * bin size in kb
let LIBZ=$LIBS/$KM#

bamCoverage --bam ${TMPDIR}/${Sample}/${Sample}.sorted.nodup.noM.bam --binSize 5 \
   --outFileFormat bigwig --smoothLength 150  \
   --scaleFactor $LIBZ \
    -o ${TMPDIR}/$Sample/$Sample.readsInPeaks.bin5.centered.smooth.150.max200.bw --numberOfProcessors ${NSLOTS} \
    --maxFragmentLength 150 \
    --centerReads --extendReads


#bamCoverage --bam ${TMPDIR}/$Sample/${Sample}.sorted.nodup.noM.bam --binSize 5 \
#    --outFileFormat bigwig --smoothLength 150  \
#    --scaleFactor $LIBZ \
#    -o ${TMPDIR}/$Sample/$Sample.readsInPeaks.bin5.centered.smooth.150.max500.bw --numberOfProcessors ${NSLOTS} \
#    --maxFragmentLength 500 \
#    --centerReads --extendReads
    #pyatac ins --bam $bam --bed $bed --

rsync -r -a -v $TMPDIR/${Sample} $FOLDERPATH/${Sample}

#mkdir -p /athena/elementolab/scratch/asd2007/Projects/DataSets/COVERAGE

#mkdir -p /athena/elementolab/scratch/asd2007/Projects/DataSets/COVERAGE/AWS.OE

#mkdir -p /athena/elementolab/scratch/asd2007/Projects/DataSets/COVERAGE/AWS.OE/${Sample}

#rsync -r -a -v $TMPDIR/${Sample}/*.bw  /athena/elementolab/scratch/asd2007/Projects/DataSets/COVERAGE/AWS.OE/${Sample}

sambamba sort --nthreads --memory-limit 30GB \
         --nthreads ${NSLOTS} --tmpdir ${TMPDIR} --out ${TMPDIR}/${Sample}/${Sample}.sorted.bam ${TMPDIR}/${Sample}/${Sample}.bam

samtools index $TMPDIR/${Sample}/${Sample}.sorted.bam 

cp $TMPDIR/${Sample}/${Sample}.sorted.bam $TMPDIR/${Sample}/${Sample}.bam

cp $TMPDIR/${Sample}/${Sample}.sorted.bam.bai $TMPDIR/${Sample}/${Sample}.bam.bai

mkdir -p ${Sample}/IDR

### IDR

source activate bds_atac_py3

~/ATACseq/getIDR.sh   $TMPDIR/${Sample}/pseudo_reps/${Sample}.nodup.pr1.tn5.pval0.1.500k.narrowPeak.gz $TMPDIR/${Sample}/pseudo_reps/${Sample}.nodup.pr2.tn5.pval0.1.500k.narrowPeak.gz ${Sample}/pseudo_reps/${Sample}.nodup.tn5.pooled.pf.pval0.1.500K.naive_overlap.narrowPeak.gz  $TMPDIR/${Sample}/IDR/${Sample}.IDR.txt


rsync -avP $TMPDIR/${Sample} $FOLDERPATH/${Sample}
rsync -r -a -v $TMPDIR/${Sample} $FOLDERPATH/${Sample}

source deactivate

#################
################
#### ATACqc ##


FINAL_BED="${PWD}/${Sample}/${Sample}.nsorted.fixmate.nodup.noM.black.Tn5.tagAlign.gz"

WORKDIR=$PWD
OUTDIR="qc"
OUTPREFIX=$Sample
INPREFIX=$Sample
#GENOME='hg19' # This is the only genome that currently works

SAMPLE="$Sample"
# Annotation FILEPATHs

#ANNOTDIR="/home/asd2007/melnick_bcell_scratch/asd2007/reference"
#ANNOTDIR="/athena/elementolab/scratch/asd2007/reference"

#DNASE_BED="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz"
#BLACKLIST_BED="${ANNOTDIR}/${GENOME}/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
#TSS_BED="${ANNOTDIR}/${GENOME}/hg19_RefSeq_stranded.bed.gz"
#REF_FASTA="${ANNOTDIR}/${GENOME}/encodeHg19Male.fa"
#PROM="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_prom_p2.bed.gz"
#ENH="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_enh_p2.bed.gz"
#REG2MAP="${ANNOTDIR}/${GENOME}/dnase_avgs_reg2map_p10_merged_named.pvals.gz"
#ROADMAP_META="${ANNOTDIR}/${GENOME}/eid_to_mnemonic.txt"
#################
################
#### ATACqc ##
# Annotation FILEPATHs
#ANNOTDIR="/home/asd2007/melnick_bcell_scratch/asd2007/reference"
#ANNOTDIR="/athena/elementolab/scratch/asd2007/reference"
#GENOME='hg19' # This is the only genome that currently works
OUTPREFIX=$Sample
INPREFIX=$Sample
SAMPLE="$Sample"
PBC="${WORKDIR}/$SAMPLE.pbc.qc"
FINAL_BAM="${Sample}/${Sample}.sorted.nodup.noM.black.bam"
FINAL_BED="${Sample}/${Sample}.nodup.tn5.tagAlign.gz"
F1="${Sample}.R1.trim.fastq.gz"
F2="${Sample}.R2.trim.fastq.gz"
ALIGNED_BAM="${Sample}/${Sample}.sorted.bam"
WORKDIR=$PWD
OUTDIR="qc"
OUTPREFIX=$Sample
INPREFIX=$Sample


export PICARD="/home/asd2007/Tools/picard/build/libs/picard.jar"
export PATH="/home/asd2007/Tools/picard/build/libs:$PATH"
#JAVA_HOME=/home/akv3001/jdk1.8.0_05
alias picard="java -Xms500m -Xmx6G -jar $PICARD"

if [ -f ${Sample}/${Sample}.picardcomplexity.qc ]; then
    echo "found picarmetrics FILEPATH"
else
    picardmetrics run -f $PICARDCONF -o ${Sample}/QCmetrics/raw ${Sample}/${Sample}.sorted.bam
    cp ${Sample}/QCmetrics/raw/*.EstimateLibraryComplexity.log ${Sample}/QCmetrics/${Sample}.picardcomplexity.qc
    cp ${Sample}/QCmetrics/raw/*.EstimateLibraryComplexity.log ${Sample}/${Sample}.picardcomplexity.qc
fi;



#python /home/asd2007/ATACseq/run_ataqc.athena.py --workdir $PWD/${Sample} \

python /home/asd2007/ATACseq/run_ataqc.athena.py --workdir $PWD/${Sample} \
    --outdir qc \
    --outprefix ${Sample} \
    --genome ${GENOME} \
    --ref ${REF_FA} --tss $TSS_ENRICH \
    --dnase ${DNASE} \
    --blacklist ${BLACK} \
    --prom $PROM \
    --enh ${ENH} \
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
   --peaks "${Sample}/peaks/${Sample}.tn5.pf.narrowPeak.gz" \
   --naive_overlap_peaks "${Sample}/pseudo_reps/${Sample}.nodup.tn5.pooled.pf.pval0.1.500K.naive_overlap.narrowPeak.gz" \
   --idr_peaks "${Sample}/IDR/${Sample}.IDR.txt.IDR0.1.filt.narrowPeak.gz"







rsync -r -a -v $TMPDIR/${Sample} $FOLDERPATH/${Sample}

rsync -r -a -v $TMPDIR/qc $FOLDERPATH/${Sample}

################
###############





if [[ NUC == 1 ]] ; then
    echo "------------------------------------------Call Nucleosomes with NucleoATAC------------------------------"
    bedtools slop -i $TMPDIR/${Sample}/${Sample}.tag.broad_peaks.broadPeak -g ${chrsz} -b 1000 > $TMPDIR/${Sample}/${Sample}.slop1k.bed
 ## save to run on pooled samples
 ## high resolution 1bp bin insertion density
   
    #pyatac ins --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.bam --out $TMPDIR/${Sample}/${Sample}.ins.smooth --smooth 151 --cores ${NSLOTS}

#
    #REF="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Homo_sapiens/UCSC/mm10/Sequence/BWAIndex/genome.fa"
    #rsync -avP /athena/elementolab/scratch/asd2007/reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.* ./
    rsync -avP $REFGen/${RG}.* ./
    nucleoatac run --bed $TMPDIR/${Sample}/${Sample}.slop1k.bed --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.bam --fasta $TMPDIR/${RG}.fa --out  $TMPDIR/${Sample}/${Sample} \
        --write_all --cores ${NSLOTS}
    #
    rsync -avP $TMPDIR/${Sample} $FOLDERPATH/${Sample}
    #rsync -r -a -v $TMPDIR/${Sample} $FOLDERPATH/${Sample}#
    ## additional tracks
    #igvtools toTDF -z 10 $TMPDIR/${Sample}/${Sample}.ins.smooth.ins.bedgraph.gz $TMPDIR/${Sample}/${Sample}.ins.smooth.ins.tdf $RG
    mv $TMPDIR/${Sample}/${Sample}.ins.smooth.ins.bedgraph.gz $TMPDIR/${Sample}/${Sample}.ins.smooth.ins.bdg.gz
    gunzip $TMPDIR/${Sample}/${Sample}.ins.smooth.ins.bdg.gz
    ~/SHELL_SCRIPTS/bdg2bw $TMPDIR/${Sample}/${Sample}.ins.smooth.ins.bdg $chrsz
    mv $TMPDIR/${Sample}/${Sample}.nucleoatac_signal.smooth.bedgraph.gz $TMPDIR/${Sample}/${Sample}.nucleoatac_signal.smooth.bdg.gz
    gunzip $TMPDIR/${Sample}/${Sample}.nucleoatac_signal.smooth.bdg.gz
    ~/SHELL_SCRIPTS/bdg2bw  $TMPDIR/${Sample}/${Sample}.nucleoatac_signal.smooth.bdg $chrsz
    rm $TMPDIR/${Sample}/${Sample}*.bdg
    ###
    rsync -avP $TMPDIR/${Sample} $FOLDERPATH/${Sample}
    rsync -r -a -v $TMPDIR/${Sample} $FOLDERPATH/${Sample}
else
    echo "---- no NucleoATAC------------------------------"
fi




#qsub  /home/asd2007/ATACseq/run_ataqc.athena.sh $GENOME



