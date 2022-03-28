#!/bin/bash

### --------------
### Whole workflow
### --------------

# 1. Build the workding directories;
# 2. Fastqc: check the quality of sequencing;
# 3. Trimgalore: filter and trim the reads and cut the adapter;
# 4. Fastqc: check the quality again;
# 4. Bowtie/Bowtie2: mapping & depended on read length;
# 5a. Sambamba: re-format sam files and remove duplicates;
# 5b. Samtools: extract uniquely-mapped reads or filter alignments with low mapping quality;
# 5c. Samtools: merge replication;
# 6. Macs2: call peak;
# 7. IDR: identify the highly reproducible peaks;
# 8. Homer/Deeptools: make vis file;

### ---------------
### Running command
### ---------------

# nohup bash ChIPseq.sh -c ChIPseq.conf &

### --------------
### Parse augments
### --------------

Help()
{
   echo -e ">-----------------------------------------------<"
   echo "Usage: $0 -c configuration file"
   echo "Example: nohup bash ChIPseq.sh -c ChIPseq.conf &"
   echo "Function: analyze ChIP-seq or Cut-Tag data: general pipeline"
   echo -e "Parameters:"
   echo -e "\t-c configuration file"
   echo -e ">-----------------------------------------------<"
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
   EGS=...
   Macs_spe=...
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
       format=`echo ${file} | awk -F'[.]' '{print $(NF-1)"."$NF}'`
       if [ "${format}" == "fastq.gz" ]; then
          prefix=`echo ${file} | sed 's,.fastq.gz,,g'`
       elif [ "${format}" == "fq.gz" ]; then
          prefix=`echo ${file} | sed 's,.fq.gz,,g'`
       fi 
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
Trimgalore(){
   indir=$1; file=$2; outdir=$3; logdir=$4; qual=$5; len=$6; ada=$7 
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   if [ ${dt} == "PE" ]; then
     for i in `seq 1 2 $(cat ${file} | wc -l)`
     do  
         read1=$(awk -v row=${i} '(NR == row){print $0}' ${file})
         r1pre=`echo ${read1} | xargs basename -s ".fq.gz"`
         read2=$(awk -v row=${i} '(NR == row+1){print $0}' ${file})
         r2pre=`echo ${read1} | xargs basename -s ".fq.gz"`
         trim_galore --paired ${indir}/${read1} ${indir}/${read2} -o ${outdir} \
                     --quality ${qual} --max_n 4 --length ${len} \
                     --${ada} --cores 16 > ${logdir}/${r1pre}_${r2pre}.log 2>&1
     done 
   elif [ ${dt} == "SE" ]; then
     cat ${file} | while read file
     do
         prefix=`echo ${file} | xargs basename -s ".fq.gz"`
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
         read1=$(awk -v row=${i} '(NR == row){print $0}' ${file})
         r1pre=`echo ${read1} | xargs basename -s ".fq.gz"`
         read2=$(awk -v row=${i} '(NR == row+1){print $0}' ${file})
         r2pre=`echo ${read1} | xargs basename -s ".fq.gz"`
         bowtie -p 16 ${GX} -1 ${indir}/${read1} -2 ${indir}/${read2} \
                -S ${outdir}/${r1pre}_${r2pre}.sam > ${logdir}/${r1pre}_${r2pre}.log 2>&1
     done
   elif [ "${dt}" == "SE" ]; then
     cat ${file} | while read file
     do
         prefix=`echo ${file} | xargs basename -s ".fq.gz"`
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
         read1=$(awk -v row=${i} '(NR == row){print $0}' ${file})
         r1pre=`echo ${read1} | xargs basename -s ".fq.gz"`
         read2=$(awk -v row=${i} '(NR == row+1){print $0}' ${file})
         r2pre=`echo ${read1} | xargs basename -s ".fq.gz"`
         bowtie2 --local -q -N 1 -p 16 --no-mixed --no-discordant -x ${GX} -1 ${indir}/${read1} -2 ${indir}/${read2} \
                 -S ${outdir}/${r1pre}_${r2pre}.sam > ${logdir}/${r1pre}_${r2pre}.log 2>&1
     done
   elif [ "${dt}" == "SE" ]; then
     cat ${file} | while read file
     do
         prefix=`echo ${file} | xargs basename -s ".fq.gz"`
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
   cat ${file} | while read bam
   do
       prefix="${bam%.*}"
       samtools view -Sh -@ 16 ${indir}/${bam} | grep -v "XS:i:" > ${outdir}/${prefix}_unique.sam
       find ${outdir} -name "*.sam" | xargs basename -s ".sam" | xargs -P 3 -I{} sambamba view -t 12 -f bam -S -o ${outdir}/{}.bam ${outdir}/{}.sam
	   find ${outdir} -name "*.bam" | xargs basename -s ".bam" | xargs -P 3 -I{} sambamba index -t 12 ${outdir}/{}.bam
       rm ${outdir}/*.sam
   done
}
MergeTwoBam(){
   indir=$1; file=$2; outdir=$3
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}
   for i in `seq 1 2 $(cat ${file} | wc -l)`
   do
       b1=$(awk -v row=${i} '(NR == row){print $0}' ${file}); b1pre="${b1%.*}"
       b2=$(awk -v row=${i} '(NR == row+1){print $0}' ${file}); b2pre="${b2%.*}"
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
   cat ${file} | while read bam
   do
       prefix="${bam%.*}"
       Factor=$(samtools idxstats ${indir}/${bam} | cut -f3 | awk -v COUNT=${num} 'BEGIN {total=0} {total += $1} END {print COUNT/total}')
       if [[ $Factor > 1 ]]; then
          echo '[ERROR]: Requested number of reads exceeds total read count in '${file}'-- exiting' && exit 1
       fi
       samtools view -s $Factor -b ${indir}/${bam} > ${outdir}/${prefix}_${num}.bam
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
# Deeptools
DeepCoverage(){
   indir=$1; file=$2; outdir=$3; logdir=$4; gs=$5
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   if [ ${dt} == "PE" ]; then
      cat ${file} | while read bam
      do
          prefix="${bam%.*}"
          bamCoverage --bam ${indir}/${bam} --outFileName ${outdir}/${prefix}.bw --outFileFormat bigwig \
                      --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize ${gs} \
                      --extendReads --ignoreDuplicates --numberOfProcessors 16 > ${logdir}/${prefix}_coverage.log 2>&1
      done
   elif [ ${dt} == "SE" ]; then
      cat ${file} | while read bam
      do
          prefix="${bam%.*}"
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
          input=$(awk -v row=${i} '(NR == row){print $0}' ${file}); input_pre="${input%.*}"
          ip=$(awk -v row=${i} '(NR == row+1){print $0}' ${file}); ip_pre="${ip%.*}"
          bamCompare -b1 ${indir}/${ip} -b2 ${indir}/${input} -o ${outdir}/${ip_pre}_vs_${input_pre}.bw --outFileFormat bigwig \
                     --operation "log2" --pseudocount 1 --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize ${gs} \
                     --scaleFactorsMethod None --extendReads --ignoreDuplicates --numberOfProcessors 16 > ${logdir}/${ip_pre}_vs_${input_pre}.log 2>&1
      done
   elif [ ${dt} == "SE" ]; then
      for i in `seq 1 2 $(cat ${file} | wc -l)`
      do
          input=$(awk -v row=${i} '(NR == row){print $0}' ${file}); input_pre="${input%.*}"
          ip=$(awk -v row=${i} '(NR == row+1){print $0}' ${file}); ip_pre="${ip%.*}"
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
# Multiqc
Multiqc(){
   indir=$1; outdir=$2; logdir=$3
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   multiqc ${indir} -o ${outdir} > ${logdir}/run.log 2>&1
}
# Homer
HomerBw(){
   indir=$1; file=$2; outdir=$3; logdir=$4; chromsize=$5
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   cat ${file} | while read bam
   do
       prefix="${bam%.*}"
       makeTagDirectory ${outdir}/${prefix}_TagDir ${indir}/${bam} > ${logdir}/${prefix}_TagDir.log 2>&1
       makeUCSCfile ${outdir}/${prefix}_TagDir -o auto -fsize 5e7 -res 1 > ${logdir}/${prefix}_bedgraph.log 2>&1
       makeUCSCfile ${outdir}/${prefix}_TagDir -o auto -bigWig ${chromsize} -fsize 1e20 > ${logdir}/${prefix}_bigwig.log 2>&1
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
Macs2WithCon(){
   indir=$1; file=$2; outdir=$3; logdir=$4; mode=$5; pvalue=$6; spe=$7; readsout=$8
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   if [ "${mode}" == "broad" ]; then
      if [ "${readsout}" == "PE" ]; then
         cat ${file} | awk '(NR >= 2){print $1}'| while read bam
         do
             input=$(cat ${file} | awk '(NR == 1){print $1}'); prefix="${bam%.*}"
             macs2 callpeak -t ${indir}/${bam} -c ${indir}/${input} -n ${prefix} \
                            -p ${pvalue} --broad --broad-cutoff ${pvalue} -g ${spe} -f BAMPE --outdir ${outdir} > ${logdir}/${prefix} 2>&1
             cut -f 1-3 ${outdir}/${prefix}_peaks.broadPeak | sort -k1,1 -k2,2n > ${outdir}/${prefix}_bpks.bed
         done
      elif [ "${readsout}" == "SE" ]; then
         cat ${file} | awk '(NR >= 2){print $1}'| while read bam
         do
             input=$(cat ${file} | awk '(NR == 1){print $1}'); prefix="${bam%.*}"
             macs2 callpeak -t ${indir}/${bam} -c ${indir}/${input} -n ${prefix} \
                            -p ${pvalue} --broad --broad-cutoff ${pvalue} -g ${spe} --nomodel --shift 37 --extsize 73 --outdir ${outdir} > ${logdir}/${prefix} 2>&1
             cut -f 1-3 ${outdir}/${prefix}_peaks.broadPeak | sort -k1,1 -k2,2n > ${outdir}/${prefix}_bpks.bed
         done
      else
         echo "Error: Unrecognized data type!"; exit 1
      fi
   elif [ "${mode}" == "narrow" ]; then
      if [ "${readsout}" == "PE" ]; then
         cat ${file} | awk '(NR >= 2){print $1}'| while read bam
         do
             input=$(cat ${file} | awk '(NR == 1){print $1}'); prefix="${bam%.*}"
             macs2 callpeak -t ${indir}/${bam} -c ${indir}/${input} -n ${prefix} \
                            -p ${pvalue} -g ${spe} -f BAMPE --outdir ${outdir} > ${logdir}/${prefix} 2>&1
             cut -f 1-3 ${outdir}/${prefix}_peaks.narrowPeak | sort -k1,1 -k2,2n > ${outdir}/${prefix}_npks.bed
         done
      elif [ "${readsout}" == "SE" ]; then
         cat ${file} | awk '(NR >= 2){print $1}'| while read bam
         do
             input=$(cat ${file} | awk '(NR == 1){print $1}'); prefix="${bam%.*}"
             macs2 callpeak -t ${indir}/${bam} -c ${indir}/${input} -n ${prefix} \
                            -p ${pvalue} -g ${spe} --nomodel --shift 37 --extsize 73 --outdir ${outdir} > ${logdir}/${prefix} 2>&1
             cut -f 1-3 ${outdir}/${prefix}_peaks.narrowPeak | sort -k1,1 -k2,2n > ${outdir}/${prefix}_npks.bed
         done
      else
         echo "Error: Unrecognized data type!"; exit 1
      fi
   else
      echo "Error: Unrecongized the peak type!"; exit 1
   fi
}
Macs2WithoutCon(){
   indir=$1; file=$2; outdir=$3; logdir=$4; mode=$5; pvalue=$6; spe=$7; readsout=$8
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   if [ "${mode}" == "broad" ]; then
      if [ "${readsout}" == "PE" ]; then
         cat ${file} | awk '{print $1}'| while read bam
         do
             prefix="${bam%.*}"
             macs2 callpeak -t ${indir}/${bam} -n ${prefix} --broad -p ${pvalue} -g ${spe} \
                            -f BAMPE --outdir ${outdir} > ${logdir}/${prefix} 2>&1
             cut -f 1-3 ${outdir}/${prefix}_peaks.broadPeak | sort -k1,1 -k2,2n > ${outdir}/${prefix}_bpks.bed
         done
      elif [ "${readsout}" == "SE" ]; then
         cat ${file} | awk '{print $1}'| while read bam
         do
             prefix="${bam%.*}"
             macs2 callpeak -t ${indir}/${bam} -n ${prefix} --broad -p ${pvalue} -g ${spe} \
                            --nomodel --shift 37 --extsize 73 --outdir ${outdir} > ${logdir}/${prefix} 2>&1
             cut -f 1-3 ${outdir}/${prefix}_peaks.broadPeak | sort -k1,1 -k2,2n > ${outdir}/${prefix}_bpks.bed
         done
      else
         echo "Error: Unrecognized data type!"; exit 1
      fi
   elif [ "${mode}" == "narrow" ]; then
      if [ "${readsout}" == "PE" ]; then
         cat ${file} | awk '{print $1}'| while read bam
         do
             prefix="${bam%.*}"
             macs2 callpeak -t ${indir}/${bam} -n ${prefix} -p ${pvalue} -g ${spe} \
                            -f BAMPE --outdir ${outdir} > ${logdir}/${prefix} 2>&1
             cut -f 1-3 ${outdir}/${prefix}_peaks.narrowPeak | sort -k1,1 -k2,2n > ${outdir}/${prefix}_npks.bed
         done
      elif [ "${readsout}" == "SE" ]; then
         cat ${file} | awk '{print $1}'| while read bam
         do
             prefix="${bam%.*}"
             macs2 callpeak -t ${indir}/${bam} -n ${prefix} -p ${pvalue} -g ${spe} \
                            --nomodel --shift 37 --extsize 73 --outdir ${outdir} > ${logdir}/${prefix} 2>&1
             cut -f 1-3 ${outdir}/${prefix}_peaks.narrowPeak | sort -k1,1 -k2,2n > ${outdir}/${prefix}_npks.bed
         done
      else
         echo "Error: Unrecognized data type!"; exit 1
      fi
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
       pk1=$(awk -v row=${i} '(NR == row){print $0}' ${file}); pk1pre="${pk1%.*}"; pk1post="${pk1##*.}"
       pk2=$(awk -v row=${i} '(NR == row+1){print $0}' ${file}); pk2pre="${pk2%.*}"; pk2post="${pk2##*.}"
       sort -k8,8nr ${indir}/${pk1} > ${outdir}/${pk1pre}_sorted.${pk1post}
       sort -k8,8nr ${indir}/${pk2} > ${outdir}/${pk2pre}_sorted.${pk2post}
       idr --samples ${outdir}/${pk1pre}_sorted.${pk1post} ${outdir}/${pk2pre}_sorted.${pk2post} \
           --input-file-type ${mode} --rank p.value --output-file ${outdir}/${pk1pre}_${pk2pre}.${mode} \
           --plot --log-output-file ${logdir}/${pk1pre}_${pk2pre}.log
       cut -f 1-3 ${outdir}/${pk1pre}_${pk2pre}.${mode} | sort -k1,1 -k2,2n > ${outdir}/${pk1pre}_${pk2pre}.bed
   done
}
# Count reads number in Bam files
CountBam(){
   indir=$1; files=$2; outdir=$3
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}
   fullname="${files##*/}"; prefix="${fullname%.*}"
   cat ${files} | while read file
   do
       if [ "${file##*.}" == "sam" ] || [ "${file##*.}" == "bam" ]; then
         total=$(samtools idxstats ${indir}/${file} | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {print total}')
         echo ${total}  >> ${outdir}/temp.txt
       else
         echo "Error: Unrecognized data type!"; exit 1
       fi
   done
   paste ${files} ${outdir}/temp.txt > ${outdir}/${prefix}_count.txt; rm ${outdir}/temp.txt
}

### ---
### Run
### ---

echo "1st step. Build working shop"
echo "Begin at $(date)"
mkdir logs results rawdata metadata scripts
cd logs && mkdir fastqc trimgalore bowtie2 samtools sambamba deeptools \
macs2 idr bedtools multiqc homer picard
cd ../results && mkdir fastqc trimgalore bowtie2 samtools sambamba deeptools \
macs2 idr bedtools multiqc homer picard && cd ..
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
Qual=20; Len=30; Ada=illumina
Trimgalore ${wd}/rawdata ${wd}/metadata/raw_fq_${dt}.txt ${wd}/results/trimgalore ${wd}/logs/trimgalore ${Qual} ${Len} ${Ada}
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
# 1. Re-format sam files and remove duplicates
SamToDedupPara ${wd}/results/bowtie2 ${wd}/results/sambamba
if [ "${dt}" == "PE" ]; then
   ls ${wd}/results/sambamba | grep "val" | grep "dedup.bam$" > ${wd}/metadata/psort_dedup_bam_${dt}.txt
else
   ls ${wd}/results/sambamba | grep "trimmed" | grep "dedup.bam$" > ${wd}/metadata/psort_dedup_bam_${dt}.txt
fi
CountBam ${wd}/results/sambamba ${wd}/metadata/psort_dedup_bam_${dt}.txt ${wd}/results/countbam
# 2. Merge bam files
MergeTwoBam ${wd}/results/sambamba ${wd}/metadata/all_Ip_dedup_bam_file.txt ${wd}/results/samtools/merge
ls ${wd}/results/samtools/merge | grep ".bam$" > ${wd}/metadata/merge_bam_file.txt
# 3. Extract unique-mapping reads from merged bam files
ExtractUnique ${wd}/results/samtools/merge ${wd}/metadata/merge_bam_file.txt \
              ${wd}/results/samtools/merge_uni_with_dis ${wd}/logs/samtools/merge_uni_with_dis
ls ${wd}/results/samtools/merge_uni_with_dis/ | grep ".bam$" > ${wd}/metadata/merged_unique_with_discordant_bam.txt
# 4. Extract reads with mapping quality larger than 10
mkdir ${wd}/results/samtools/mapq10
find ${wd}/results/samtools/merge/ -name "*.bam" | xargs basename -s ".bam" | \
     xargs -P 3 -I{} samtools view -h -q 10 -b -o ${wd}/results/samtools/mapq10/{}_mapq10.bam ${wd}/results/samtools/merge/{}.bam
find ${wd}/results/samtools/mapq10 -name "*mapq10.bam" | xargs basename -s ".bam" | xargs -P 3 -I{} sambamba index -t 12 ${wd}/results/samtools/mapq10/{}.bam
echo "Finish at $(date)"



echo "7th step. Visulization"
echo "Begin at $(date)"
# 1. make bw files with deeptools: best-mapping
DeepCoverage ${wd}/results/sambamba ${wd}/metadata/merge_bam_file.txt \
             ${wd}/results/deeptools/bigwig/merged_bam ${wd}/logs/deeptools/bigwig/merged_bam ${EGS}
# 2. make bw files with deeptools: unique-mapping reads
DeepCoverage ${wd}/results/samtools/merge_uni_with_dis ${wd}/metadata/merged_unique_with_discordant_bam.txt \
             ${wd}/results/deeptools/bigwig/merge_uni_with_dis ${wd}/logs/deeptools/bigwig/merge_uni_with_dis ${EGS}
# 3. make bw files with deeptools: reads with mapping quality larger than 10
ls ${wd}/results/samtools/mapq10 | grep ".bam$"  > ${wd}/metadata/mapq10_bam.txt
FunDeepCoverage ${wd}/results/samtools/mapq10/ ${wd}/metadata/mapq10_bam.txt \
                ${wd}/results/deeptools/bigwig/mapq10_bam ${wd}/logs/deeptools/bigwig/mapq10_bam 2913022398
echo "Finish at $(date)"



echo "8th step. Call peaks"
echo "Begin at $(date)"
for sample in h4C93 h8C93 hICM93 hM93 hTE93
do
   Macs2WithCon ${wd}/results/sambamba ${wd}/metadata/${sample}_callpeak.txt \
                ${wd}/results/macs2/p_0.05 ${wd}/logs/macs2/p_0.05 broad 0.05 ${Macs_spe} ${dt}
done
echo "Finish at $(date)"



echo "9th step. Identify reproducible peaks"
echo "Begin at $(date)"
PeakIDR ${wd}/results/macs2/p_0.05 ${wd}/metadata/peak_for_idr.txt ${wd}/results/idr/p_0.05 ${wd}/logs/idr/p_0.05 broadPeak
echo "Finish at $(date)"
