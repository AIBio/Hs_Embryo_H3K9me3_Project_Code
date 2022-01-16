#!/bin/bash

### --------------
### Whole workflow
### --------------

# 1. Build the workding directories (required);
# 2. Fastqc: check the quality of sequencing (required);
# 3. Trimgalore: filter and trim the reads and cut the adapter (required);
# 4. Fastqc: check the quality again (required);
# 5. Bowtie/Bowtie2: mapping & depended on read length (required);
# 6a. Sambamba: re-format sam files and remove duplicates (required);
# 6b. Samtools: extract uniquely-mapped reads (optional & strict processing method);
# 6c. Samtools: merge replication or subsample alignments in fix number (optional);
# 7. Macs2: call peak (required, ATAC-seq !!!);
# 8. IDR: identify the highly reproducible peaks (required if biological replicates exist);
# 9. Homer/Deeptools/Picard: make vis file for qc ... (optional);

### ---------------
### Running command
### ---------------

# nohup bash ATACv1.sh -c ATACv1.conf &

### --------------
### Parse augments
### --------------

Help()
{
   echo -e ">--------------------------------------------------<"
   echo "Usage: $0 -c configuration file"
   echo "Example: nohup bash ATACv1.sh -c ATACv1.conf &"
   echo "Function: analyze ATAC-seq data with normal alignment"
   echo -e "Parameters:"
   echo -e "\t-c configuration file"
   echo -e ">--------------------------------------------------<"
   exit 1
}
while getopts "c:h" opt
do
   case "$opt" in
      c ) cf="$OPTARG" ;;
      h ) Help ;;
   esac
done
if [ -z "${cf}" ]; then
   echo "Some or all of the parameters are empty"
   Help
fi

echo "> ---------------------------- <"
# 1. workding directory
wd=`grep -v "^#" $cf | awk '{if(NR==1) print $0}'`;         echo -e "[1] working directory: ${wd}"
# 2. full path of raw data
rawdata=`grep -v "^#" $cf | awk '{if(NR==2) print $0}'`;    echo -e "[2] directory of raw fastq files: ${rawdata}"
# 3. txt file that records sample names
samples=`grep -v "^#" $cf | awk '{if(NR==3) print $0}'`;    echo -e "[3] sample file: ${samples}"
# 4. data type: support PE and SE data
dt=`grep -v "^#" $cf | awk '{if(NR==4) print $0}'`;         echo -e "[4] sample file: ${dt}"
# 5. species: only support human, mouse and rat
sp=`grep -v "^#" $cf | awk '{if(NR==5) print $0}'`;         echo -e "[5] species: ${sp}"
# 6. index of mapping: Bowtie2 index
index=`grep -v "^#" $cf | awk '{if(NR==6) print $0}'`;      echo -e "[6] genome index: ${index}"
# 7. chromosome size
cs=`grep -v "^#" $cf | awk '{if(NR==7) print $0}'`;         echo -e "[7] chromosome size: ${cs}"
# add new parameters
if [ "${sp}" == "human" ]; then
   EGS=2913022398
   Macs_spe=hs
elif [ "${sp}" == "mouse" ]; then
   EGS=2652783500
   Macs_spe=mm
elif [ "${sp}" == "rat" ]; then
   EGS=2729860805
   Macs_spe=2729860805
else
   echo "Unrecognized or unsupported species"; exit 1
fi
echo "> ---------------------------- <"

### --------------------
### Define the functions
### --------------------

# Create soft links for raw data
DataLink(){
   indir=$1; sample=$2; outdir=$3
   [ ! -d ${outdir} ] && mkdir -p ${outdir}
   ls ${indir} | grep -f ${sample} - | while read file
   do
       prefix="${file%%.*}"
       echo $prefix
       ln -s ${indir}/${file} ${outdir}/${prefix}.fq.gz
   done
}
# Rename files in batch
Rename(){
   indir=$1; Key=$2; Str1=$3; Str2=$4
   for file in ${indir}/*.${Key}; do mv "$file" "${file/${Str1}/${Str2}}"; done
}
# Fastqc
QcReads()
{
   indir=$1; outdir=$2; logdir=$3
   [ ! -d ${outdir} ] && mkdir -p ${outdir}; [ ! -d ${logdir} ] && mkdir -p ${logdir}
   find ${indir} -name "*.fq.gz" | xargs -P 10 -I{} fastqc {} -o ${outdir} > ${logdir}/run.log 2>&1
}
# Multiqc
ParseQc()
{
   indir=$1; outdir=$2; logdir=$3
   [ ! -d ${outdir} ] && mkdir -p ${outdir}; [ ! -d ${logdir} ] && mkdir -p ${logdir}
   multiqc ${indir} -o ${outdir} > ${logdir}/run.log 2>&1
}
# Trimgalore
TrimReads(){
   indir=$1; file=$2; outdir=$3; logdir=$4; qual=$5; len=$6; ada=$7
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   if [ ${dt} == "PE" ]; then
     for i in `seq 1 2 $(cat ${file} | wc -l)`
     do
         read1=$(awk -v row=${i} '(NR == row){print $0}' ${file}); r1pre="${read1%%.*}"
         read2=$(awk -v row=${i} '(NR == row+1){print $0}' ${file}); r2pre="${read2%%.*}"
         trim_galore --paired ${indir}/${read1} ${indir}/${read2} -o ${outdir} \
                     --quality ${qual} --max_n 4 --length ${len} \
                     --${ada} --cores 16 > ${logdir}/${r1pre}_${r2pre}.log 2>&1
     done
   elif [ ${dt} == "SE" ]; then
     cat ${file} | while read file
     do
         prefix="${file%%.*}"
         trim_galore ${indir}/${file} -o ${outdir} --quality ${qual} --max_n 4 \
                     --length ${len} --${ada} --cores 16 > ${logdir}/${prefix}.log 2>&1
     done
   else
     echo "Unrecognized data type!"; exit 1
   fi
}
# Bowtie
BowtieAlign(){
   indir=$1; file=$2; outdir=$3; logdir=$4; GX=$5
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   if [ "${dt}" == "PE" ]; then
     for i in `seq 1 2 $(cat ${file} | wc -l)`
     do
         read1=$(awk -v row=${i} '(NR == row){print $0}' ${file}); r1pre="${read1%%.*}"
         read2=$(awk -v row=${i} '(NR == row+1){print $0}' ${file}); r2pre="${read2%%.*}"
         bowtie -p 16 ${GX} -1 ${indir}/${read1} -2 ${indir}/${read2} \
                -S ${outdir}/${r1pre}_${r2pre}.sam > ${logdir}/${r1pre}_${r2pre}.log 2>&1
     done
   elif [ "${dt}" == "SE" ]; then
     cat ${file} | while read file
     do
         prefix="${file%%.*}"
         bowtie -p 16 ${GX} ${indir}/${file} \
                -S ${outdir}/${prefix}.sam > ${logdir}/${prefix}.log 2>&1
     done
   else
     echo "Unrecognized data type!"; exit 1
   fi
}
# Bowtie2 (customize the parameters by actual situation: Mismatch = 1 & no-mixed & no-discordant ==> by default)
Bowtie2Align(){
   indir=$1; file=$2; outdir=$3; logdir=$4; GX=$5
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   if [ "${dt}" == "PE" ]; then
     for i in `seq 1 2 $(cat ${file} | wc -l)`
     do
         read1=$(awk -v row=${i} '(NR == row){print $0}' ${file}); r1pre="${read1%%.*}"
         read2=$(awk -v row=${i} '(NR == row+1){print $0}' ${file}); r2pre="${read2%%.*}"
         bowtie2 --local -q -N 1 -p 16 --no-mixed --no-discordant -x ${GX} -1 ${indir}/${read1} -2 ${indir}/${read2} \
                 -S ${outdir}/${r1pre}_${r2pre}.sam > ${logdir}/${r1pre}_${r2pre}.log 2>&1
     done
   elif [ "${dt}" == "SE" ]; then
     cat ${file} | while read file
     do
         prefix="${file%%.*}"
         bowtie2 --local -q -N 1 -p 16 -x ${GX} -U ${indir}/${file} \
                 -S ${outdir}/${prefix}.sam > ${logdir}/${prefix}.log 2>&1
     done
   else
     echo "Error: Unrecognized data type!"; exit 1
   fi
}
# Samtools
ExtractUnique(){
   indir=$1; file=$2; outdir=$3; logdir=$4
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   if [ ${dt} == "PE" ]; then
      cat ${file} | while read bam
      do
          prefix="${bam%%.*}"
          samtools view -Sh -f 2 -F 256 -@ 16 ${indir}/${bam} | grep -v "XS:i:" > ${outdir}/${prefix}_unique.sam
          python /home/yhw/document/ScriptsLibrary/Python/CheckSamPairs.py < ${outdir}/${prefix}_unique.sam > ${outdir}/${prefix}_uniquePE.sam
          find ${outdir} -type f ! -name '*uniquePE.sam' -delete
      done
   elif [ ${dt} == "SE" ]; then
      cat ${file} | while read bam
      do
           prefix="${bam%%.*}"
           samtools view -Sh -F 256 -@ 16 ${indir}/${bam} | grep -v "XS:i:" > ${outdir}/${prefix}_uniqueSE.sam
       done
   else
      echo "Error: Unrecognized data type!"; exit 1
   fi
}
MergeTwoBam(){
   indir=$1; file=$2; outdir=$3
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}
   for i in `seq 1 2 $(cat ${file} | wc -l)`
   do
       b1=$(awk -v row=${i} '(NR == row){print $0}' ${file}); b1pre="${b1%%.*}"
       b2=$(awk -v row=${i} '(NR == row+1){print $0}' ${file}); b2pre="${b2%%.*}"
       samtools merge -@ 16 -O BAM ${outdir}/${b1pre}_${b2pre}_merged.bam ${indir}/${b1} ${indir}/${b2}
   done
   find ${outdir} -name "*merged.bam" | xargs -P 6 -I{} sambamba index -t 12 {}
}
MergeMoreBam(){
   intxt=$1; outdir=$2; prefix=$3
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}
   samtools merge -@ 10 -b ${intxt} ${outdir}/${prefix}_merge.bam
   sambamba index -t 12 ${outdir}/${prefix}_merge.bam
}
SubsetBam(){
   indir=$1; file=$2; outdir=$3; num=$4
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}
   cat ${file} | while read file
   do
       prefix="${file%%.*}"
       Factor=$(samtools idxstats ${indir}/${file} | cut -f3 | awk -v COUNT=${num} 'BEGIN {total=0} {total += $1} END {print COUNT/total}')
       if [[ $Factor > 1 ]]; then
          echo '[ERROR]: Requested number of reads exceeds total read count in '${file}'-- exiting' && exit 1
       fi
       samtools view -s $Factor -b ${indir}/${file} > ${outdir}/${prefix}_${num}.bam
       sambamba index -t 12 ${outdir}/${prefix}_${num}.bam
   done
}
# Sambamba
SamToDedupPara(){
   indir=$1; outdir=$2
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   find ${indir} -name "*.sam" | xargs basename -s ".sam" | xargs -P 3 -I{} sambamba view -t 12 -f bam -S -o ${outdir}/{}.bam ${indir}/{}.sam
   find ${outdir} -name "*.bam" | xargs basename -s ".bam" | xargs -P 3 -I{} sambamba sort -o ${outdir}/{}_psort.bam ${outdir}/{}.bam
   find ${outdir} -name "*psort.bam" | xargs basename -s ".bam" | xargs -P 3 -I{} sambamba markdup -t 12 -r ${outdir}/{}.bam ${outdir}/{}_dedup.bam
   find ${outdir} -name "*dedup.bam" | xargs basename -s ".bam" | xargs -P 3 -I{} sambamba index -t 12 ${outdir}/{}.bam
   find ${outdir} -type f ! -name '*psort*' -delete
}
# Picard
InsertSize(){
   indir=$1; file=$2; outdir=$3; logdir=$4
   [ ! -d "${outdir}" ] && mkdir ${outdir}; [ ! -d "${logdir}" ] && mkdir ${logdir}
   cat ${file} | while read bam
   do
       prefix="${bam%%.*}"
       picard CollectInsertSizeMetrics I=${indir}/${bam} O=${outdir}/${prefix}_insert_size_matrix.txt \
                                       H=${outdir}/${prefix}_insert_size_histogram.pdf M=0.5 > ${logdir}/${prefix}.log 2>&1
   done
}
# Deeptools
DeepCoverage(){
   indir=$1; file=$2; outdir=$3; logdir=$4; gs=$5
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   if [ ${dt} == "PE" ]; then
      cat ${file} | while read bam
      do
          prefix="${bam%%.*}"
          bamCoverage --bam ${indir}/${bam} --outFileName ${outdir}/${prefix}.bw --outFileFormat bigwig \
                      --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize ${gs} \
                      --extendReads --ignoreDuplicates --numberOfProcessors 16 > ${logdir}/${prefix}_coverage.log 2>&1
      done
   elif [ ${dt} == "SE" ]; then
      cat ${file} | while read bam
      do
          prefix="${bam%%.*}"
          bamCoverage --bam ${indir}/${bam} --outFileName ${outdir}/${prefix}.bw --outFileFormat bigwig \
                      --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize ${gs} \
                      --ignoreDuplicates --numberOfProcessors 16 > ${logdir}/${prefix}_coverage.log 2>&1
      done
   else
      echo "Error: Unrecognized data type!"; exit 1
   fi
}
BamCompare(){
   indir=$1; file=$2; outdir=$3; logdir=$4; gs=$5
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   if [ ${dt} == "PE" ]; then
      for i in `seq 1 2 $(cat ${file} | wc -l)`
      do
          input=$(awk -v row=${i} '(NR == row){print $0}' ${file}); input_pre="${input%%.*}"
          ip=$(awk -v row=${i} '(NR == row+1){print $0}' ${file}); ip_pre="${ip%%.*}"
          bamCompare -b1 ${indir}/${ip} -b2 ${indir}/${input} -o ${outdir}/${ip_pre}_vs_${input_pre}.bw --outFileFormat bigwig \
                     --operation "log2" --pseudocount 1 --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize ${gs} \
                     --scaleFactorsMethod None --extendReads --ignoreDuplicates --numberOfProcessors 16 > ${logdir}/${ip_pre}_vs_${input_pre}.log 2>&1
      done
   elif [ ${dt} == "SE" ]; then
      for i in `seq 1 2 $(cat ${file} | wc -l)`
      do
          input=$(awk -v row=${i} '(NR == row){print $0}' ${file}); input_pre="${input%%.*}"
          ip=$(awk -v row=${i} '(NR == row+1){print $0}' ${file}); ip_pre="${ip%%.*}"
          bamCompare -b1 ${indir}/${ip} -b2 ${indir}/${input} -o ${outdir}/${ip_pre}_vs_${input_pre}.bw --outFileFormat bigwig \
                     --operation "log2" --pseudocount 1 --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize ${gs} \
                     --scaleFactorsMethod None --ignoreDuplicates --numberOfProcessors 16 > ${logdir}/${ip_pre}_vs_${input_pre}.log 2>&1
      done
   else
      echo "Error: Unrecognized data type!"; exit 1
   fi
}
GlobalCoverageNor(){
   indir=$1; outdir=$2; logdir=$3; key=$4; prefix=$5
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   bwlist=$(ls ${indir}/*${key}*.bw | tr '\n' ' ')
   multiBigwigSummary bins -bs 10000 -n 0 -p 16 -b ${bwlist} -o ${outdir}/${prefix}_readNor.npz \
                           --outRawCounts ${outdir}/${prefix}_readNor.tab > ${logdir}/${prefix}_nor.log 2>&1
}
RegionCoverageNor(){
   indir=$1; outdir=$2; logdir=$3; key=$4; region=$5
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   bwlist=$(ls ${indir}/*${key}*.bw | tr '\n' ' ')
   fileName="${region##*/}"; prefix="${fileName%.*}"
   multiBigwigSummary BED-file --BED ${region} -p 16 -b ${bwlist} -o ${outdir}/${prefix}_readNor.npz \
                               --outRawCounts ${outdir}/${prefix}_readNor.tab > ${logdir}/${prefix}_nor.log 2>&1
}
# Homer
HomerBw(){
   indir=$1; file=$2; outdir=$3; logdir=$4
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   cat ${file} | while read bam
   do
       prefix="${bam%.*}"
       makeTagDirectory ${outdir}/${prefix}_TagDir ${indir}/${bam} > ${logdir}/${prefix}_TagDir.log 2>&1
       makeUCSCfile ${outdir}/${prefix}_TagDir -o auto -fsize 5e7 -res 1 > ${logdir}/${prefix}_bedgraph.log 2>&1
       makeUCSCfile ${outdir}/${prefix}_TagDir -o auto -bigWig ${cs} -fsize 1e20 > ${logdir}/${prefix}_bigwig.log 2>&1
   done
}
HomerMotif(){
   indir=$1; file=$2; outdir=$3; logdir=$4; gv=$5
   cat ${file} | while read region
   do
       prefix="${region%.*}"
       findMotifsGenome.pl ${indir}/${region} ${gv} ${outdir} > ${logdir}/${prefix}.log 2>&1
   done
}
HomerPeakAnno(){
   indir=$1; file=$2; outdir=$3; gv=$4
   cat ${file} | while read region
   do
       prefix="${region%.*}"
       annotatePeaks.pl ${indir}/${region} ${gv} > ${outdir}/${prefix}_homer_annotation.txt
   done
}
# Macs2
CallPeak(){
   indir=$1; file=$2; outdir=$3; logdir=$4; mode=$5; pvalue=$6; spe=$7
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   if [ "${mode}" == "broad" ]; then
      cat ${file} | awk '{print $1}'| while read bam
      do
          prefix="${bam%%.*}"
          macs2 callpeak -t ${indir}/${bam} -n ${prefix} --broad --broad-cutoff ${pvalue} -g ${spe} \
                         -f BAMPE --outdir ${outdir} > ${logdir}/${prefix} 2>&1
          cut -f 1-3 ${outdir}/${prefix}_peaks.broadPeak | sort -k1,1 -k2,2n > ${outdir}/${prefix}_bpks.bed
      done
   elif [ "${mode}" == "narrow" ]; then
      cat ${file} | awk '{print $1}'| while read bam
      do
          prefix="${bam%%.*}"
          macs2 callpeak -t ${indir}/${bam} -n ${prefix} -p ${pvalue} -g ${spe} \
                         -f BAMPE --outdir ${outdir} > ${logdir}/${prefix} 2>&1
          cut -f 1-3 ${outdir}/${prefix}_peaks.narrowPeak | sort -k1,1 -k2,2n > ${outdir}/${prefix}_npks.bed
      done
   else
      echo "Error: Unrecongized the peak type !"; exit 1
   fi
}
# IDR
PeakIDR(){
   indir=$1; file=$2; outdir=$3; logdir=$4; mode=$5
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   for i in `seq 1 2 $(cat ${file} | wc -l)`
   do
       pk1=$(awk -v row=${i} '(NR == row){print $0}' ${file}); pk1pre="${pk1%%.*}"; pk1post="${pk1##*.}"
       pk2=$(awk -v row=${i} '(NR == row+1){print $0}' ${file}); pk2pre="${pk2%%.*}"; pk2post="${pk2##*.}"
       sort -k8,8nr ${indir}/${pk1} > ${outdir}/${pk1pre}_sorted.${pk1post}
       sort -k8,8nr ${indir}/${pk2} > ${outdir}/${pk2pre}_sorted.${pk2post}
       idr --samples ${outdir}/${pk1pre}_sorted.${pk1post} ${outdir}/${pk2pre}_sorted.${pk2post} \
           --input-file-type ${mode} --rank p.value --output-file ${outdir}/${pk1pre}_${pk2pre}.${mode} \
           --plot --log-output-file ${logdir}/${pk1pre}_${pk2pre}.log
       cut -f 1-3 ${outdir}/${pk1pre}_${pk2pre}.${mode} | sort -k1,1 -k2,2n > ${outdir}/${pk1pre}_${pk2pre}.bed
   done
}

### ---
### Run
### ---

echo "1st step. Build working shop"
echo "Begin at $(date)"
mkdir logs results rawdata metadata scripts
cd logs && mkdir fastqc trimgalore bowtie2 samtools sambamba deeptools \
macs2 idr bedtools multiqc ngsplot snpsplit homer picard
cd ../results && mkdir fastqc trimgalore bowtie2 samtools sambamba deeptools \
macs2 idr bedtools multiqc ngsplot snpsplit homer picard && cd ..
DataLink ${rawdata} ${samples} ${wd}/rawdata
echo "Finish at $(date)"



echo "2nd step. Check the quality of raw data"
echo "Begin at $(date)"
QcReads ${wd}/rawdata ${wd}/results/fastqc/rawdata ${wd}/logs/fastqc/rawdata
ParseQc ${wd}/results/fastqc/rawdata ${wd}/results/multiqc/rawdata ${wd}/logs/multiqc RawData
if [ "${dt}" == "PE" ]; then
   ls ${wd}/rawdata | grep "_[1-2]" > ${wd}/metadata/raw_fq_${dt}.txt
else
   ls ${wd}/rawdata | grep -v "_[1-2]" > ${wd}/metadata/raw_fq_${dt}.txt
fi
echo "Finish at $(date)"



echo "3rd step. Filter and trim the reads"
echo "Begin at $(date)"
Qual=20; Len=30; Ada=nextera
TrimReads ${wd}/rawdata ${wd}/metadata/raw_fq_${dt}.txt ${wd}/results/trimgalore ${wd}/logs/trimgalore ${Qual} ${Len} ${Ada}
if [ "${dt}" == "PE" ]; then
   ls ${wd}/results/trimgalore | grep "fq.gz$" | grep "val" > ${wd}/metadata/trimgalore_fq_${dt}.txt
else
   ls ${wd}/results/trimgalore | grep "fq.gz$" | grep "trimmed" > ${wd}/metadata/trimgalore_fq_${dt}.txt
fi
echo "Finish at $(date)"



echo "4th step. Check the quality of raw data again"
echo "Begin at $(date)"
QcReads ${wd}/results/trimgalore ${wd}/results/fastqc/trimgalore ${wd}/logs/fastqc/trimgalore
ParseQc ${wd}/results/fastqc/trimgalore ${wd}/results/multiqc/trimgalore ${wd}/logs/multiqc TrimData
echo "Finish at $(date)"



echo "5th step. Align read to genome"
echo "Begin at $(date)"
Bowtie2Align ${wd}/results/trimgalore ${wd}/metadata/trimgalore_fq_${dt}.txt ${wd}/results/bowtie2 ${wd}/logs/bowtie2 ${index}
ParseQc ${wd}/logs/bowtie2 ${wd}/results/multiqc/bowtie2 ${wd}/logs/multiqc Bowtie2Alignment
if [ "${dt}" == "PE" ]; then
   ls ${wd}/results/bowtie2 | grep "val" | grep "sam$" > ${wd}/metadata/sam_${dt}.txt
else
   ls ${wd}/results/bowtie2 | grep "trimmed" | grep "sam$" > ${wd}/metadata/sam_${dt}.txt
fi
echo "Finish at $(date)"



echo "6th step. Process the sam files"
echo "Begin at $(date)"
### >>> Extract unique-mapping reads (optional)
#ExtractUnique ${wd}/results/bowtie2 ${wd}/metadata/sam_${dt}.txt ${wd}/results/samtools ${wd}/logs/samtools ${dt}
#if [ "${dt}" == "PE" ]; then
#   ls ${wd}/results/samtools | grep "val" | grep "unique" > ${wd}/metadata/unique_sam_${dt}.txt
#else
#   ls ${wd}/results/samtools | grep "trimmed" | grep "unique" > ${wd}/metadata/unique_sam_${dt}.txt
#fi
### >>> Re-format sam files and remove duplicates
SamToDedupPara ${wd}/results/bowtie2 ${wd}/results/sambamba
if [ "${dt}" == "PE" ]; then
   ls ${wd}/results/sambamba | grep "val" | grep "dedup.bam$" > ${wd}/metadata/psort_dedup_bam_${dt}.txt
else
   ls ${wd}/results/sambamba | grep "trimmed" | grep "dedup.bam$" > ${wd}/metadata/psort_dedup_bam_${dt}.txt
fi
echo "Finish at $(date)"



echo "7th step. Visulization"
echo "Begin at $(date)"

### >>> make bw files with homer
HomerBw ${wd}/results/sambamba ${wd}/metadata/psort_dedup_bam_${dt}.txt ${wd}/results/homer ${wd}/logs/homer
### >>> make bw files with deeptools: effective genome size (GRCh38-2913022398) (GRCm38-2652783500)
DeepCoverage ${wd}/results/sambamba ${wd}/metadata/psort_dedup_bam_${dt}.txt \
             ${wd}/results/deeptools/coverage ${wd}/logs/deeptools/coverage ${EGS}
### >>> collect insert size
InsertSize ${wd}/results/sambamba ${wd}/metadata/psort_dedup_bam_${dt}.txt ${wd}/results/picard/insert_size ${wd}/logs/picard/insert_size
echo "Finish at $(date)"



echo "8th step. Call peaks"
echo "Begin at $(date)"
CallPeak ${wd}/results/sambamba ${wd}/metadata/psort_dedup_bam_${dt}.txt \
         ${wd}/results/macs2/broad_p0.05 ${wd}/logs/macs2/broad_p0.05 broad 0.05 ${Macs_spe}
CallPeak ${wd}/results/sambamba ${wd}/metadata/psort_dedup_bam_${dt}.txt \
         ${wd}/results/macs2/narrow_p0.05 ${wd}/logs/macs2/narrow_p0.05 narrow 0.05 ${Macs_spe}
echo "Finish at $(date)"



#echo "9th step. Identify reproducible peaks"
#echo "Begin at $(date)"
#ls ${wd}/results/macs2/broad_p0.05 | grep "broadPeak$" > ${wd}/metadata/macs2_bpks.txt
#PeakIDR ${wd}/results/macs2/broad_p0.05 ${wd}/metadata/macs2_bpks.txt ${wd}/results/idr/broad_p0.05 ${wd}/logs/idr/broad_p0.05 broadPeak
#ls ${wd}/results/macs2/narrow_p0.05 | grep "narrowPeak$" > ${wd}/metadata/macs2_npks.txt
#PeakIDR ${wd}/results/macs2/narrow_p0.05 ${wd}/metadata/macs2_npks.txt ${wd}/results/idr/narrow_p0.05 ${wd}/logs/idr/narrow_p0.05 narrowPeak
#echo "Finish at $(date)"
