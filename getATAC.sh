#!/bin/bash -l
#$ -N ATACseq
#$ -j y
#$ -m a
#$ -cwd
#$ -l zenodotus=true
#$ -l os=rhel6.3
#$ -M ashley.doane@gmail.com
#$ -l h_rt=23:00:00
#$ -pe smp 8-16
#$ -l h_vmem=12G
#$ -R y


########### SETTINGS ###########
TRIM=1 # 
NUC=0 # run nucleoatac, extended runtime required
BT2ALN=1 # bt2 recommnded!
##############################


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


rsync -r -v -a -z $path/$file/* ./

rsync -r -v -a -z $path/$file/*.sra ./
#rsync -r -v -a -z $path/$file/ ./
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
    #BLACK="/home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/Anshul_Hg19UltraHighSignalArtifactRegions.bed"
    #chrsz="/home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19.genome.chrom.sizes"
    fetchChromSizes hg19 > hg19.chrom.sizes
    chrsz= $PWD/hg19.chrom.sizes
    RG="hg19"
    SPEC="hs"
    REFGen="/home/asd2007/melnick_bcell_scratch/asd2007/bin/bcbio/genomes/Hsapiens/hg19/seq/"
    rsync -av /home/asd2007/Scripts/picardmetrics.conf ./

elif [ $gtf_path == 2 ]
then
	gtf="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Mus_UCSC_ref.gtf"
	REF="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex"
    REFbt2="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index"
    BLACK="/home/asd2007/dat02/asd2007/Reference/mm10-blacklist.bed"
    RG="mm10"
    SPEC="mm"
    REFGen="/home/asd2007/melnick_bcell_scratch/asd2007/bin/bcbio/genomes/Mmusculus/mm10/seq/"
    #chrsz="/home/asd2007/melnick_bcell_scratch/asd2007/Reference/mm10.genome.chrom.sizes"
    fetchChromSizes mm10 > mm10.chrom.sizes
    chrsz= $PWD/mm10.chrom.sizes
    rsync -av /home/asd2007/Scripts/picardmetrics.Mouse.conf ./
else
	gtf="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/mm10_UCSC_ref.gtf"
	REF="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Homo_sapiens/UCSC/mm10/Sequence/BWAIndex/genome.fa"
fi

mkdir ref

rsync -avP ${REF} $TMPDIR/

rsync -avP ${REFbt2} $TMPDIR/

mkdir ${TMPDIR}/${Sample}

echo "ls of pwd is"
ls -lrth



if [ $TRIM == 1 ]
then
    echo "Trimming adapter sequences, with command..."

#    cat *_val_1.fq.gz > ${Sample}.R1.fastq.gz
 #   cat *_val_2.fq.gz > ${Sample}.R2.fastq.gz
    cat *_R2* >  ${Sample}.R2.fastq.gz
    cat *_R1* > ${Sample}.R1.fastq.gz
    #cat *R2* >  ${Sample}.R2.fastq.gz
    trimAdapters.py -a ${Sample}.R1.fastq.gz -b ${Sample}.R2.fastq.gz
    echo "completed trimming"
    pigz -p $NSLOTS -c $TMPDIR/${Sample}.R1.trim.fastq > $TMPDIR/${Sample}.R1.trim.fastq.gz
    pigz -p $NSLOTS -c $TMPDIR/${Sample}.R2.trim.fastq > $TMPDIR/${Sample}.R2.trim.fastq.gz
    rsync -av $TMPDIR/${Sample}.R*.trim.fastq.gz ${path}/${Sample}
else
    echo "will not perform adapter sequence trimming of reads"
    #gunzip *.gz
      #mv ${Sample}.R1.fastq $TMPDIR/${Sample}.R1.trim.fastq
      #mv ${Sample}.R2.fastq $TMPDIR/${Sample}.R2.trim.fastq
    # mv $TMPDIR/R2/${Sample}.R2.fastq $TMPDIR/R2/${Sample}.R2.trim.fastq
    # mv $TMPDIR/R1/${Sample}.R1.fastq $TMPDIR/R1/${Sample}.R1.trim.fastq
fi


#bowtie2  -X 2000 -p ${NSLOTS} -x $TMPDIR/Bowtie2Index/genome -1 $TMPDIR/R1/${Sample}.R1.fq -2 $TMPDIR/R2/${Sample}.R2.fq -S $TMPDIR/${Sample}/${Sample}.bt2.sam



if [ $BT2ALN == 1 ]
then
    echo "----------bowtie2 aligning-------------"
    bowtie2 -X 2000 --threads ${NSLOTS} -x $TMPDIR/Bowtie2Index/genome  \
            -1 $TMPDIR/${Sample}.R1.trim.fastq.gz -2 $TMPDIR/${Sample}.R2.trim.fastq.gz 2> ${Sample}/${Sample}.align.log| \
        samtools view -bS - > $TMPDIR/${Sample}/${Sample}.bam
   # samtools sort  $TMPDIR/${Sample}.bam -o $TMPDIR/${Sample}/${Sample}.sorted.bam
    #cp $TMPDIR/${Sample}/${Sample}.sorted.bam  $TMPDIR/${Sample}/${Sample}.bt2.sorted.bam``
else
    echo "----------bwa-mem aligning-------------"
   # echo "aligning : $TMPDIR/${Sample}.R1.trim.fq ,  $TMPDIR/${Sample}.R2.trim.fq using bwa-mem.."
   # bwa mem -t ${NSLOTS} -M $TMPDIR/BWAIndex/genome.fa $TMPDIR/${Sample}.R1.trim.fastq $TMPDIR/${Sample}.R2.trim.fastq | samtools view -bS - >  $TMPDIR/${Sample}.bam
fi

echo "----------Processing alignment and filtering for duplicates and mitochondrial mapping reads-----------------"

processBamAlignment.sh $TMPDIR/${Sample}/${Sample}.bam

samtools sort -n $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam -o $TMPDIR/${Sample}/${Sample}.nsorted.nodup.noM.black.bam
FINAL_BAM='${TMPDIR}/${Sample}/${Sample}.sorted.nodup.noM.black.bam'

echo "--------------------------Tn5 adjusted bedfile for MACS2 peak calling---------------------------"

samtools fixmate $TMPDIR/${Sample}/${Sample}.nsorted.nodup.noM.black.bam $TMPDIR/${Sample}/${Sample}.nsorted.fixmate.nodup.noM.black.bam

convertBAMtoBED.sh $TMPDIR/${Sample}/${Sample}.nsorted.fixmate.nodup.noM.black.bam

samtools view -H $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.bam | grep chr | grep -v chrM | /home/ole2001/PERL_SCRIPTS/columns.pl 1 2 | sed 's/SN://' | sed 's/LN://' > chrom.sizes

#TMPDIR/${Sample}/${Sample}.nsorted.fixmate.nodup.noM.black.bam


echo "------------------------------------------Call Peaks with MACS2 ON PSEUDOREPLICTES--------------------------------------------"

callPeaks.sh $TMPDIR/${Sample}/${Sample}.nsorted.fixmate.nodup.noM.black.bedpe.gz ${Sample} hs

#adjustedBed="/home/ole2001/PROGRAMS/SOFT/bedtools2/bin/slopBed -i $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bed -g sizes -l 75 -r -75 -s"

#macs2 callpeak -t <(${adjustedBed}) -f BED -n $TMPDIR/${Sample}/${Sample}.broad -g 2.7e9 -p 1e-3 --nomodel --shift 75 -B --SPMR --broad --keep-dup all


echo "------------------------------------------Call Peaks with MACS2--------------------------------------------"

macs2 callpeak -t  $TMPDIR/${Sample}/${Sample}.nsorted.fixmate.nodup.noM.black.Tn5.tagAlign.gz -f BED -n $TMPDIR/${Sample}/${Sample}.tag.narrow -g $SPEC  --nomodel --shift -75 --extsize 150 --keep-dup all --bdg --SPMR --call-summits -p 1e-3

sort -k 8gr,8gr $TMPDIR/${Sample}/${Sample}.tag.narrow_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | gzip -nc > $TMPDIR/${Sample}/${Sample}.tn5.pf.narrowPeak.gz

macs2 callpeak -t  $TMPDIR/${Sample}/${Sample}.nsorted.fixmate.nodup.noM.black.Tn5.tagAlign.gz -f BED -n $TMPDIR/${Sample}/${Sample}.tag.broad -g $SPEC  --nomodel --shift -75 --extsize 150 --keep-dup all --broad --broad-cutoff 0.1

sort -k 8gr,8gr $TMPDIR/${Sample}/${Sample}.tag.broad_peaks.broadPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | gzip -nc > $TMPDIR/${Sample}/${Sample}.tn5.pf.broadPeak.gz

sort -k 14gr,14gr $TMPDIR/${Sample}/${Sample}.tag.broad_peaks.gappedPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | gzip -nc > $TMPDIR/${Sample}/${Sample}.tn5.pf.gappedPeak.gz


#macs2 callpeak -t  $TMPDIR/${Sample}/${Sample}.nsorted.fixmate.nodup.noM.black.bedpe.gz -f BEDPE -n $TMPDIR/${Sample}/${Sample}.bedpe.narrow -g $SPEC  --nomodel --shift -75 --extsize 150 --keep-dup all --call-summits -p 1e-3

#macs2 callpeak -t $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam -f BAMPE -n $TMPDIR/${Sample}/${Sample}.narrow -g $SPEC --nomodel --shift -75 --extsize 150 --keep-dup all --call-summits -p 1e-3

#macs2 callpeak -t $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam -f BAMPE -n $TMPDIR/${Sample}/${Sample}.broad -g $SPEC  --nomodel --shift -75 --extsize 150 --keep-dup all --broad --broad-cutoff 0.1


#macs2 callpeak -t ${TMPDIR}/${Sample}/${Sample}.sorted.nodup.noM.frag123.bam -f BAMPE -n $TMPDIR/${Sample}/${Sample}.frag123.narrow -g hs --nomodel --shift -75 --extsize 150 --keep-dup all --call-summits -p 1e-3

## lenient threshold for IDR
#macs2 callpeak -t $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam -f BAMPE -n $TMPDIR/${Sample}/${Sample}.narrow.p0.01 -g hs --nomodel --shift -75 --extsize 150 --keep-dup all --call-summits -p 1e-1

#macs2 callpeak -t $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam -f BAMPE -n $TMPDIR/${Sample}/${Sample}.broad.p0.1 -g hs  --nomodel --shift -100 --extsize 200 --keep-dup all --broad --broad-cutoff 0.1
#rsync -avP /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.* ./
#rsync -avP /home/asd2007/dat02/asd2007/Reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.* ./
prefix=$TMPDIR/${Sample}/${Sample}.narrow

#macs2 bdgcmp -t ${prefix}_treat_pileup.bdg -c ${prefix}_control_lambda.bdg
#    --o-prefix ${prefix} -m FE
#slopBed -i ${prefix}_FE.bdg -g $chrsz -b 0 > ${prefix}_clip.bdg

#bedClip ${prefix}_clip.bdg $chrsz {$prefix}_fc.bdg
#rm -f ${prefix}_FE.bdg
#rm -f ${prefix}_clip.bdg

#sort -k 1,1 -k 2,2n ${prefix}_fc.bdg > ${prefix}_fc_srt.bdg
#bedGraphToBigWig ${prefix}_fc_srt.bdg $chrsz ${prefix}_fc.bigwig
#rm -f ${prefix}_fc_srt.bdg ${prefix}_fc.bdg


#echo "--- deeptools coverage -rpkm --"

#bamCoverage --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam --binSize 5 \
#    --outFileFormat bigwig --smoothLength 200  \
#    --normalizeUsingRPKM \
#    -o $TMPDIR/${Sample}/${Sample}.bamCov.rpkm.bw --centerReads --numberOfProcessors ${NSLOTS}



echo "----------- compute QC stats  ------------------"

mkdir -p ${Sample}/QCmetrics

mkdir -p ${Sample}/QCmetrics/raw

mkdir -p ${Sample}/QCmetrics/filtered



pyatac sizes --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam --upper 1000 --out $TMPDIR/${Sample}/QCmetrics/${Sample}.filtered


picardmetrics run -o $TMPDIR/${Sample}/QCmetrics/raw $TMPDIR/${Sample}/${Sample}.sorted.bam

picardmetrics run -o $TMPDIR/${Sample}/QCmetrics/filtered $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam

cp $TMPDIR/${Sample}/QCmetrics/filtered/*.EstimateLibraryComplexity.log $TMPDIR/${Sample}/QCmetrics/${Sample}.picardcomplexity.qc

cp $TMPDIR/${Sample}/QCmetrics/filtered/*.EstimateLibraryComplexity.log $TMPDIR/${Sample}/${Sample}.picardcomplexity.qc

samtools flagstat  $TMPDIR/${Sample}/${Sample}.sorted.bam > $TMPDIR/${Sample}/QCmetrics/raw/${Sample}.sorted.flagstat.txt

samtools flagstat $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.bam > $TMPDIR/${Sample}/QCmetrics/filtered/${Sample}.sorted.nodup.noM.flagstat.txt

getchrM $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam > $TMPDIR/${Sample/}QCmetrics/${Sample}.sorted.nodup.noM.black.bam.chrM.txt


getchrM $TMPDIR/${Sample}/${Sample}.sorted.bam > $TMPDIR/$Sample}/QCmetrics/${Sample}.sorted.bam.chrM.txt


echo "----------- compute fragment midpoint coverage per bp ------------------"
#pyatac cov --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam --out $TMPDIR/${Sample}/${Sample}.pyatac.cov --cores ${NSLOTS}

#pyatac cov --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam --out $TMPDIR/${Sample}/${Sample}.pyatac.125bp.cov --cores ${NSLOTS} --upper 125

bamCoverage --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam --binSize 5 \
    --outFileFormat bigwig --smoothLength 150 \
    --normalizeUsingRPKM \
    --maxFragmentLength 300 \
    -o $TMPDIR/${Sample}/${Sample}.smooth151.center.extend.fpkm.max300.bw --centerReads --extendReads --numberOfProcessors ${NSLOTS}


bamCoverage --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam --binSize 5 \
            --outFileFormat bigwig --smoothLength 150 \
            --normalizeUsingRPKM \
            --maxFragmentLength 150 \
            -o $TMPDIR/${Sample}/${Sample}.smooth151.center.extend.fpkm.max115.bw --centerReads --extendReads --numberOfProcessors ${NSLOTS}


#bamCoverage --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam --binSize 20 \
#    --outFileFormat bigwig --smoothLength 150 \
#    --normalizeUsingRPKM \
#    -o $TMPDIR/${Sample}/${Sample}.bin20.smooth150.fpkm.bw --numberOfProcessors ${NSLOTS}

rsync -avP $TMPDIR/${Sample} $path/${Sample}
rsync -r -a -v $TMPDIR/${Sample} $path/${Sample}





BED="/home/asd2007/dat02/asd2007/Projects/DataSets/atacData/ATAC-AML/AML/cellAtlas.bed"
#BED="/home/asd2007/dat02/asd2007/Projects/DataSets/atacData/atacCN1/project/cellAtlas.CN1.bed"


pyatac counts --bam  $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam --bed ${BED} --out $TMPDIR/${Sample}/${Sample}.ins --upper 200
##bedtools multicov -split -bams $BAMS -bed $BED

LIBS=$(zcat ${Sample}/${Sample}.ins.counts.txt.gz  | awk '{sum +=$1} END  {printf sum}')

#let KM=2000000 #1 million * bin size in kb
#let LIBZ=$LIBS/$KM



#bamCoverage --bam  $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam --binSize 20 \
#    --outFileFormat bigwig --smoothLength 150  \
#    --scaleFactor $LIBZ \
#    -o ${TMPDIR}/${Sample}/$Sample.readsInPeaks.bin20.centered.smooth.150.max150f.bw --numberOfProcessors ${NSLOTS} \
#    --maxFragmentLength 150 \
#    --centerReads --extendReads



#bamCoverage --bam ${TMPDIR}/${Sample}/${Sample}.sorted.nodup.noM.black.bam --binSize 20 \
#    --outFileFormat bigwig --smoothLength 150  \
#    --scaleFactor $LIBZ \
#    -o ${TMPDIR}/$Sample/$Sample.readsInPeaks.bin20.centered.smooth.150.max500f.bw --numberOfProcessors ${NSLOTS} \
#    --maxFragmentLength 500 \
#    --centerReads --extendReads


let KM=500000 #1 million * bin size in kb
let LIBZ=$LIBS/$KM

bamCoverage --bam ${TMPDIR}/${Sample}/${Sample}.sorted.nodup.noM.black.bam --binSize 5 \
    --outFileFormat bigwig --smoothLength 150  \
    --scaleFactor $LIBZ \
    -o ${TMPDIR}/$Sample/$Sample.readsInPeaks.bin5.centered.smooth.150.max200.bw --numberOfProcessors ${NSLOTS} \
    --maxFragmentLength 200 \
    --centerReads --extendReads


#bamCoverage --bam ${TMPDIR}/$Sample/${Sample}.sorted.nodup.noM.black.bam --binSize 5 \
#    --outFileFormat bigwig --smoothLength 150  \
#    --scaleFactor $LIBZ \
#    -o ${TMPDIR}/$Sample/$Sample.readsInPeaks.bin5.centered.smooth.150.max500.bw --numberOfProcessors ${NSLOTS} \
#    --maxFragmentLength 500 \
#    --centerReads --extendReads
    #pyatac ins --bam $bam --bed $bed --

rsync -avP $TMPDIR/${Sample} $path/${Sample}
rsync -r -a -v $TMPDIR/${Sample} $path/${Sample}

mkdir -p /home/asd2007/melnick_bcell_scratch/asd2007/COVERAGE/AWS.OE/${Sample}
rsync -r -a -v $TMPDIR/${Sample}/*.bw  /home/asd2007/melnick_bcell_scratch/asd2007/COVERAGE/AWS.OE/${Sample}

#path="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Jon_bwa_mm10_output"

#rsync -r -v $TMPDIR/${Sample}* /zenodotus/dat01/melnick_bcell_scratch/asd2007/COVERAGE/ATAC/test/peaksNorm/

















Sample="sample_LY7"
FINAL_BAM="${Sample}/${Sample}.sorted.nodup.noM.black.bam"
FINAL_BED="${Sample}/${Sample}.nsorted.fixmate.nodup.noM.black.Tn5.tagAlign.gz"

WORKDIR="${Sample}"
OUTDIR="${WORKDIR}/QCmetrics/atacqc"
OUTPREFIX="${Sample}"
INPREFIX="$Sample"
GENOME='hg19' # This is the only genome that currently works

SAMPLE="$Sample"
# Annotation files
ANNOTDIR="/home/asd2007/melnick_bcell_scratch/asd2007/Reference"
X="/home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/"
DNASE_BED="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz"
BLACKLIST_BED="/home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
TSS_BED="${ANNOTDIR}/${GENOME}/hg19_RefSeq_stranded.bed.gz"
REF_FASTA="${ANNOTDIR}/${GENOME}/encodeHg19Male.fa"
PROM="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_prom_p2.bed.gz"
ENH="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_enh_p2.bed.gz"
REG2MAP="${ANNOTDIR}/${GENOME}/dnase_avgs_reg2map_p10_merged_named.pvals.gz"
ROADMAP_META="${ANNOTDIR}/${GENOME}/eid_to_mnemonic.txt"
PBC="${WORKDIR}/$SAMPLE.pbc.qc"
F1="${Sample}.R1.trim.fastq.gz"
F2="${Sample}.R2.trim.fastq.gz"


python $HOME/dat02/asd2007/TOOLS/ataqc/run_ataqc.py --workdir ${Sample} \
    --outdir ${Sample} \
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
    --fastq1 ${Sample}.R1.trim.fastq.gz \
    --fastq2 ${Sample}.R2.trim.fastq.gz \
    --alignedbam ${Sample}/${Sample}.bam \
    --alignmentlog ${Sample}/${Sample}.align.log \
    --coordsortbam ${Sample}/${Sample}.bam \
    --duplog ${Sample}/${Sample}.sorted.dup.qc \
    --pbc ${Sample}/${Sample}.sorted.pbc.qc \
    --finalbam ${FINAL_BAM} \
    --finalbed ${FINAL_BED} \
    --bigwig $Sample/$Sample.readsInPeaks.bin5.centered.smooth.150.max200.bw \
    --peaks ${Sample}/${Sample}.tn5.pf.narrowPeak.gz \
    --naive_overlap_peaks ${Sample}/${Sample}.tn5_pooled.pf.pval0.1.500K.naive_overlap.narrowPeak.gz


















#
#
#
#
#
#
#
#
#
#
#
#
#
#
#FINAL_BAM="${TMPDIR}/${Sample}/${Sample}.sorted.nodup.noM.black.bam"
#FINAL_BED="$TMPDIR/${Sample}/${Sample}.nsorted.fixmate.nodup.noM.black.Tn5.tagAlign.gz"
#
#WORKDIR="${TMPDIR}/${Sample}"
#OUTDIR="${WORKDIR}/QCmetrics/atacqc"
#OUTPREFIX="${Sample}"
#INPREFIX="$Sample"
#GENOME='hg19' # This is the only genome that currently works
#
#SAMPLE="$Sample"
## Annotation files
#ANNOTDIR="/home/asd2007/melnick_bcell_scratch/asd2007/Reference"
#X="/home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/"
#DNASE_BED="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz"
#BLACKLIST_BED="/home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
#TSS_BED="${ANNOTDIR}/${GENOME}/hg19_RefSeq_stranded.bed.gz"
#REF_FASTA="${ANNOTDIR}/${GENOME}/encodeHg19Male.fa"
#PROM="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_prom_p2.bed.gz"
#ENH="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_enh_p2.bed.gz"
#REG2MAP="${ANNOTDIR}/${GENOME}/dnase_avgs_reg2map_p10_merged_named.pvals.gz"
#ROADMAP_META="${ANNOTDIR}/${GENOME}/eid_to_mnemonic.txt"
#PBC="${WORKDIR}/$SAMPLE.pbc.qc"
#F1="$TMPDIR/${Sample}.R1.trim.fastq.gz"
#F2="$TMPDIR/${Sample}.R2.trim.fastq.gz"
#
#
#
#
#
#
#
#
#
#python $HOME/dat02/asd2007/TOOLS/ataqc/run_ataqc.py --workdir ${TMPDIR}/${Sample} \
#    --outdir ${TMPDIR}/${Sample} \
#    --outprefix ${Sample} \
#    --genome hg19 \
#    --ref /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/encodeHg19Male/encodeHg19Male.fa \
#    --tss /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/hg19_RefSeq_stranded.bed.gz \
#    --dnase /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz \
#    --blacklist /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz \
#    --prom /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/reg2map_honeybadger2_dnase_prom_p2.bed.gz \
#    --enh /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/reg2map_honeybadger2_dnase_enh_p2.bed.gz \
#    --reg2map /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/dnase_avgs_reg2map_p10_merged_named.pvals.gz \
#    --meta /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/eid_to_mnemonic.txt \
#    --fastq1 $TMPDIR/${Sample}.R1.trim.fastq.gz \
#    --fastq2 $TMPDIR/${Sample}.R2.trim.fastq.gz \
#    --alignedbam ${TMPDIR}/${Sample}/${Sample}.bam \
#    --aligmmentlog ${TMPDIR}/${Sample}/${Sample}.align.log \
#    --coordsortbam ${TMPDIR}/${Sample}/${Sample}.bam \
#    --duplog ${TMPDIR}/${Sample}/${Sample}.sorted.dup.qc \
#    --pbc ${TMPDIR}/${Sample}/${Sample}.sorted.pbc.qc \
#    --finalbam ${FINAL_BAM} \
#    --finalbed ${FINAL_BED} \
#    --bigwig ${TMPDIR}/$Sample/$Sample.readsInPeaks.bin5.centered.smooth.150.max200.bw \
#    --peaks $TMPDIR/${Sample}/${Sample}.tn5.pf.narrowPeak.gz \
#    --naive_overlap_peaks $TMPDIR/${Sample}/${Sample}.tn5_pooled.pf.pval0.1.500K.naive_overlap.narrowPeak.gz 
    #--idr_peaks \



                   # --naive_overlap_peaks /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/macs2/overlap/Sample_GCB_138.R1.trim.PE2SE.nodup.tn5_pooled.pf.pval0.1.500K.naive_overlap.narrowPeak.gz \
		                 # --idr_peaks /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/idr/optimal_set/ppr.IDR0.1.filt.narrowPeak.gz \





if [ NUC == 1 ]
then
    echo "------------------------------------------Call Nucleosomes with NucleoATAC------------------------------"
    bedtools slop -i $TMPDIR/${Sample}/${Sample}.broad.broadPeak -g ${chrsz} -b 1000 > $TMPDIR/${Sample}/${Sample}.slop1k.bed
 ## save to run on pooled samples
 ## high resolution 1bp bin insertion density
   
    #pyatac ins --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam --out $TMPDIR/${Sample}/${Sample}.ins.smooth --smooth 151 --cores ${NSLOTS}

#
    #REF="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Homo_sapiens/UCSC/mm10/Sequence/BWAIndex/genome.fa"
    #rsync -avP /home/asd2007/dat02/asd2007/Reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.* ./
    rsync -avP $REFGen/${RG}.* ./
    nucleoatac run --bed $TMPDIR/${Sample}/${Sample}.slop1k.bed --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam --fasta $TMPDIR/${RG}.fa --out  $TMPDIR/${Sample}/${Sample} \
        --write_all --cores ${NSLOTS}
    #
    rsync -avP $TMPDIR/${Sample} $path/${Sample}
    #rsync -r -a -v $TMPDIR/${Sample} $path/${Sample}#
    ## additional tracks
    #igvtools toTDF -z 10 $TMPDIR/${Sample}/${Sample}.ins.smooth.ins.bedgraph.gz $TMPDIR/${Sample}/${Sample}.ins.smooth.ins.tdf $RG
    mv $TMPDIR/${Sample}/${Sample}.ins.smooth.ins.bedgraph.gz $TMPDIR/${Sample}/${Sample}.ins.smooth.ins.bdg.gz
    gunzip $TMPDIR/${Sample}/${Sample}.ins.smooth.ins.bdg.gz
    ~/SHELL_SCRIPTS/bdg2bw $TMPDIR/${Sample}/${Sample}.ins.smooth.ins.bdg $chrsz
    mv $TMPDIR/${Sample}/${Sample}.nucleoatac_signal.smooth.bedgraph.gz $TMPDIR/${Sample}/${Sample}.nucleoatac_signal.smooth.bdg.gz
    gunzip $TMPDIR/${Sample}/${Sample}.nucleoatac_signal.smooth.bdg.gz
    ~/SHELL_SCRIPTS/bdg2bw  $TMPDIR/${Sample}/${Sample}.nucleoatac_signal.smooth.bdg $chrsz
    #rm $TMPDIR/${Sample}/${Sample}*.bdg
    ###
    rsync -avP $TMPDIR/${Sample} $path/${Sample}
    rsync -r -a -v $TMPDIR/${Sample} $path/${Sample}
else
    echo "---- no NucleoATAC------------------------------"

    rsync -avP $TMPDIR/${Sample} $path/${Sample}
    rsync -r -a -v $TMPDIR/${Sample} $path/${Sample}
fi








fetchChromSizes hg19 > hg19.chrom.sizes

prefix=$TMPDIR/${Sample}/${Sample}.tag.narrow
#macs2 bdgcmp -t ${prefix}_treat_pileup.bdg -c ${prefix}_control_lambda.bdg -o ${prefix}_FE.bdg -m FE

macs2 bdgcmp -t ${prefix}_treat_pileup.bdg -c ${prefix}_control_lambda.bdg \
    --o-prefix ${prefix} -m FE

slopBed -i ${prefix}_FE.bdg -g $chrsz -b 0 | bedClip stdin hg19.chrom.sizes ${prefix}_fc.bedgraph
rm -f ${prefix}_FE.bdg

sort -k1,1 -k2,2n ${prefix}_fc.bedgraph > ${prefix}_fc.srt.bedgraph

bedGraphToBigWig  ${prefix}_fc.srt.bedgraph $chrsz ${prefix}_fc.bw





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


rsync -avP $TMPDIR/${Sample}* $path/${Sample}/
rsync -r -a -v $TMPDIR/${Sample} $path/${Sample}/






##
##		sys macs2 bdgcmp -t "$prefix"_treat_pileup.bdg -c "$prefix"_control_lambda.bdg --outdir $peak_o_dir -o "$prefix_basename"_FE.bdg -m FE
##
##		//# Remove coordinates outside chromosome sizes (stupid MACS2 bug)
##		sys slopBed -i "$prefix"_FE.bdg -g $chrsz -b 0 |   awk '{if ($3 != -1) print $0}' |  bedClip stdin $chrsz $fc_bedgraph
##		sys rm -f "$prefix"_FE.bdg
##		
##		//# Convert bedgraph to bigwig
##		sys sort -k1,1 -k2,2n $fc_bedgraph > $fc_bedgraph_srt
##		sys bedGraphToBigWig $fc_bedgraph_srt $chrsz $fc_bigwig
##		sys rm -f $fc_bedgraph $fc_bedgraph_srt
##
##		//===========================================
##		//# For -log10(p-value) signal tracks
##		//============================================
##		
##		sys chipReads=$(zcat $tag | wc -l | awk '{printf "%f", $1/1000000}')
##		//# Compute sval = min(no. of reads in ChIP, no. of reads in control) / 1,000,000
##		sys $sval_line
##
##		sys macs2 bdgcmp -t "$prefix"_treat_pileup.bdg -c "$prefix"_control_lambda.bdg --outdir $peak_o_dir -o "$prefix_basename"_ppois.bdg -m ppois -S "${sval}"
##
##		//# Remove coordinates outside chromosome sizes (stupid MACS2 bug)
##		sys slopBed -i "$prefix"_ppois.bdg -g $chrsz -b 0 |   awk '{if ($3 != -1) print $0}' |  bedClip stdin $chrsz $peak_o_dir/$prefix_basename.pval.signal.bedgraph
##		sys rm -rf "$prefix"_ppois.bdg
##		
##		//# Convert bedgraph to bigwig
##		sys sort -k1,1 -k2,2n $pval_bedgraph > $pval_bedgraph_srt
##		sys bedGraphToBigWig $pval_bedgraph_srt $chrsz $pval_bigwig
##		sys rm -f $pval_bedgraph $pval_bedgraph_srt
##
##		sys rm -f "$prefix"_treat_pileup.bdg "$prefix"_control_lambda.bdg








