#!/bin/bash

### ---------------
### Whole procedure
### ---------------

# 1. Build the working directory (required);
# 2. Fastqc + MultiQC: check the quality of sequencing (required);
# 3. Trimgalore + Fastqc + MultiQC: trim reads and QC (required);
# 4a. STAR: align the reads to genome by STAR (required);
# 4b. STAR: align the reads to genome to quantify repeat in subfamily level (optional);
# 5. Sambamba: re-format sam files and remove duplicate (optional);
# 6. Deeptools & Homer: make the visulization file of UCSC (required);
# 7a. FeatureCounts: quantify gene and repeat locus expression levels (required);
# 7b. TeToolKit: quantify the repeat element in name level (optional);

### ---------------
### Running command
### ---------------

# nohup bash RNAseqv1.sh -c RNAseqv1.conf &

### --------------
### Parse augments
### --------------
Help()
{
   echo -e ">--------------------------------------------------------<"
   echo "Usage: $0 -c configuration file"
   echo "Example: nohup bash RNAseqv1.sh -c RNAseqv1.conf &"
   echo "Function: analyze RNA-seq data to quantify gene and repeat."
   echo -e "Parameters:"
   echo -e "\t-c configuration file"
   echo -e ">--------------------------------------------------------<"
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
wd=`grep -v "^#" $cf | awk '{if(NR==1) print $0}'`;            echo -e "[1] working directory: ${wd}"
# 2. full path of raw data
rawdata=`grep -v "^#" $cf | awk '{if(NR==2) print $0}'`;       echo -e "[2] directory of raw fastq files: ${rawdata}"
# 3. txt file that records sample names
samples=`grep -v "^#" $cf | awk '{if(NR==3) print $0}'`;       echo -e "[3] sample file: ${samples}"
# 4. data type: support PE and SE data
dt=`grep -v "^#" $cf | awk '{if(NR==4) print $0}'`;            echo -e "[4] data type: ${dt}"
# 5. species: only support human, mouse
sp=`grep -v "^#" $cf | awk '{if(NR==5) print $0}'`;            echo -e "[5] species: ${sp}"
# 6. index of mapping for gene
star_gene=`grep -v "^#" $cf | awk '{if(NR==6) print $0}'`;     echo -e "[6] STAR genome index for gene is ${star_gene}"
# 7. index of mapping for repeat
star_repeat=`grep -v "^#" $cf | awk '{if(NR==7) print $0}'`;   echo -e "[7] STAR genome index for repeat is ${star_repeat}"
# 8. gene GTF file
gene_gtf=`grep -v "^#" $cf | awk '{if(NR==8) print $0}'`;      echo -e "[8] gene GTF annotation file: ${gene_gtf}"
# 9. repeat SAF file
repeat_saf=`grep -v "^#" $cf | awk '{if(NR==9) print $0}'`;    echo -e "[9] repeat SAF annotation file: ${repeat_saf}"
# 10. repeat GTF file
repeat_gtf=`grep -v "^#" $cf | awk '{if(NR==10) print $0}'`;   echo -e "[10] repeat GTF annotation file: ${repeat_gtf}"
# 11. adapter type: illumina/nextera
adapter=`grep -v "^#" $cf | awk '{if(NR==11) print $0}'`;      echo -e "[11] adapter type: ${adapter}"
echo "> ---------------------------- <"

### ---------
### Functions
### ---------

# Create soft links for raw data
DataLink()
{
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
TrimReads()
{
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
     while read file
     do
         prefix=`echo ${file} | xargs basename -s ".fq.gz"`
         trim_galore ${indir}/${file} -o ${outdir} --quality ${qual} --max_n 4 \
                     --length ${len} --${ada} --cores 16 > ${logdir}/${prefix}.log 2>&1
     done < ${file}
   else
     echo "Unrecognized data type!"; exit 1
   fi
}
# Star
StarGene()
{
   indir=$1; file=$2; outdir=$3; logdir=$4; indexdir=$5
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   if [ ${dt} == "PE" ]; then
     for i in `seq 1 2 $(cat ${file} | wc -l)`
     do
         r1=$(awk -v row=${i} '(NR == row){print $0}' ${file})
         r1name=`echo ${r1} | xargs basename -s ".fq.gz"`
         r2=$(awk -v row=${i} '(NR == row+1){print $0}' ${file})
         r2name=`echo ${r2} | xargs basename -s ".fq.gz"`
         STAR --runMode alignReads --runThreadN 16 --readFilesCommand zcat --outBAMsortingThreadN 10 \
              --outSAMattributes NH HI NM MD XS nM AS MD jM jI MC ch \
              --outFileNamePrefix ${outdir}/${r1name}_${r2name} --outSAMtype BAM SortedByCoordinate \
              --outSAMmultNmax -1 --outFilterMultimapNmax 10 --genomeDir ${indexdir} \
              --readFilesIn ${indir}/${r1},${indir}/${r2} > ${logdir}/${r1name}_${r2name}.log 2>&1
     done
   elif [ ${dt} == "SE" ]; then
     while read fq
     do
         prefix=`echo ${fq} | xargs basename -s ".fq.gz"`
         STAR --runMode alignReads --runThreadN 16 --readFilesCommand zcat --outBAMsortingThreadN 10 \
              --outSAMattributes NH HI NM MD XS nM AS MD jM jI MC ch \
              --outFileNamePrefix ${outdir}/${prefix} --outSAMtype BAM SortedByCoordinate \
              --outSAMmultNmax -1 --outFilterMultimapNmax 10 --genomeDir ${indexdir} \
              --readFilesIn ${indir}/${fq} > ${logdir}/${prefix}.log 2>&1
     done < ${file}
   else
     echo "Unrecognized data type!"; exit 1
   fi
}
StarRandomMus()
{
   indir=$1; file=$2; outdir=$3; logdir=$4; indexdir=$5
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   if [ ${dt} == "PE" ]; then
     for i in `seq 1 2 $(cat ${file} | wc -l)`
     do
         r1=$(awk -v row=${i} '(NR == row){print $0}' ${file})
         r1name=`echo ${r1} | xargs basename -s ".fq.gz"`
         r2=$(awk -v row=${i} '(NR == row+1){print $0}' ${file})
         r2name=`echo ${r2} | xargs basename -s ".fq.gz"`
         STAR --runMode alignReads --runThreadN 16 --readFilesCommand zcat --outBAMsortingThreadN 12 \
              --outSAMattributes NH HI NM MD XS nM AS MD jM jI MC ch \
              --outFilterMultimapNmax 5000 --outSAMmultNmax 1 --outFilterMismatchNmax 3 --outMultimapperOrder Random --winAnchorMultimapNmax 5000 \
              --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 350 --seedSearchStartLmax 30 --alignTranscriptsPerReadNmax 30000 \
              --alignWindowsPerReadNmax 30000  --alignTranscriptsPerWindowNmax 300 --seedPerReadNmax 3000 --seedPerWindowNmax 300 --seedNoneLociPerWindow 1000 \
              --outFileNamePrefix ${outdir}/${r1name}_${r2name} --outSAMtype BAM SortedByCoordinate \
              --genomeDir ${indexdir} --readFilesIn ${indir}/${r1},${indir}/${r2} > ${logdir}/${r1name}_${r2name}.log 2>&1
     done
   elif [ ${dt} == "SE" ]; then
     while read fq
     do
         prefix=`echo ${fq} | xargs basename -s ".fq.gz"`
         STAR --runMode alignReads --runThreadN 16 --readFilesCommand zcat --outBAMsortingThreadN 12 \
              --outSAMattributes NH HI NM MD XS nM AS MD jM jI MC ch \
              --outFilterMultimapNmax 5000 --outSAMmultNmax 1 --outFilterMismatchNmax 3 --outMultimapperOrder Random --winAnchorMultimapNmax 5000 \
              --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 350 --seedSearchStartLmax 30 --alignTranscriptsPerReadNmax 30000 \
              --alignWindowsPerReadNmax 30000  --alignTranscriptsPerWindowNmax 300 --seedPerReadNmax 3000 --seedPerWindowNmax 300 --seedNoneLociPerWindow 1000 \
              --outFileNamePrefix ${outdir}/${prefix} --outSAMtype BAM SortedByCoordinate \
              --genomeDir ${indexdir} --readFilesIn ${indir}/${fq} > ${logdir}/${prefix}.log 2>&1
     done < ${file}
   else
     echo "Unrecognized data type!"; exit 1
   fi
}
StarRandomHs()
{
   indir=$1; file=$2; outdir=$3; logdir=$4; indexdir=$5
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   if [ ${dt} == "PE" ]; then
     for i in `seq 1 2 $(cat ${file} | wc -l)`
     do  
         r1=$(awk -v row=${i} '(NR == row){print $0}' ${file})
         r1name=`echo ${r1} | xargs basename -s ".fq.gz"`
         r2=$(awk -v row=${i} '(NR == row+1){print $0}' ${file})
         r2name=`echo ${r2} | xargs basename -s ".fq.gz"`
         STAR --runMode alignReads --runThreadN 16 --readFilesCommand zcat --outBAMsortingThreadN 12 \
              --outSAMattributes NH HI NM MD XS nM AS MD jM jI MC ch \
              --outFilterMultimapNmax 1000 --outSAMmultNmax 1 --outFilterMismatchNmax 3 --outMultimapperOrder Random \
              --winAnchorMultimapNmax 1000 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 350 \
              --outFileNamePrefix ${outdir}/${r1name}_${r2name} --outSAMtype BAM SortedByCoordinate \
              --genomeDir ${indexdir} --readFilesIn ${indir}/${r1},${indir}/${r2} > ${logdir}/${r1name}_${r2name}.log 2>&1
     done 
   elif [ ${dt} == "SE" ]; then
     while read fq
     do  
         prefix=`echo ${fq} | xargs basename -s ".fq.gz"`
         STAR --runMode alignReads --runThreadN 16 --readFilesCommand zcat --outBAMsortingThreadN 12 \
              --outSAMattributes NH HI NM MD XS nM AS MD jM jI MC ch \
              --outFilterMultimapNmax 1000 --outSAMmultNmax 1 --outFilterMismatchNmax 3 --outMultimapperOrder Random \
              --winAnchorMultimapNmax 1000 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 350 \
              --outFileNamePrefix ${outdir}/${prefix} --outSAMtype BAM SortedByCoordinate \
              --genomeDir ${indexdir} --readFilesIn ${indir}/${fq} > ${logdir}/${prefix}.log 2>&1
     done < ${file}
   else
     echo "Unrecognized data type!"; exit 1
   fi
}
StarMulti()
{
   indir=$1; file=$2; outdir=$3; logdir=$4; indexdir=$5
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   if [ ${dt} == "PE" ]; then
     for i in `seq 1 2 $(cat ${file} | wc -l)`
     do
         r1=$(awk -v row=${i} '(NR == row){print $0}' ${file})
         r1name=`echo ${r1} | xargs basename -s ".fq.gz"`
         r2=$(awk -v row=${i} '(NR == row+1){print $0}' ${file})
         r2name=`echo ${r2} | xargs basename -s ".fq.gz"`
         STAR --runMode alignReads --runThreadN 16 --readFilesCommand zcat --outBAMsortingThreadN 10 \
              --outSAMattributes NH HI NM MD XS nM AS MD jM jI MC ch \
              --outFileNamePrefix ${outdir}/${r1name}_${r2name} --outSAMtype BAM SortedByCoordinate \
              --winAnchorMultimapNmax 200 --outFilterMultimapNmax 100 --genomeDir ${indexdir} \
              --readFilesIn ${indir}/${r1},${indir}/${r2} > ${logdir}/${r1name}_${r2name}.log 2>&1
     done
   elif [ ${dt} == "SE" ]; then
     while read fq
     do
         prefix=`echo ${fq} | xargs basename -s ".fq.gz"`
         STAR --runMode alignReads --runThreadN 16 --readFilesCommand zcat --outBAMsortingThreadN 10 \
              --outSAMattributes NH HI NM MD XS nM AS MD jM jI MC ch \
              --outFileNamePrefix ${outdir}/${prefix} --outSAMtype BAM SortedByCoordinate \
              --winAnchorMultimapNmax 200 --outFilterMultimapNmax 100 --genomeDir ${indexdir} \
              --readFilesIn ${indir}/${file} > ${logdir}/${prefix}.log 2>&1
     done < ${file}
   else
     echo "Unrecognized data type!"; exit 1
   fi
}
# Sambamba
BamToDedupPara()
{
   indir=$1; outdir=$2; logdir=$3
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   find ${indir} -name "*psorted.bam" | xargs basename -s ".bam" | parallel -P 3 --gun "sambamba markdup -r -t 12 ${indir}/{}.bam ${outdir}/{}_dedup.bam > ${logdir}/{}_dedup.log 2>&1"
   find ${outdir} -name "*dedup.bam" | xargs basename -s ".bam" | parallel -P 3 --gun "sambamba index -t 12 ${outdir}/{}.bam > ${logdir}/{}_index.log 2>&1"
}
# Subread
CountGene()
{
   indir=$1; gtf=$2; outdir=$3; logdir=$4
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   bamlist=$(ls ${indir}/*.bam | tr '\n' ' '); echo ${bamlist}
   featureCounts -a ${gtf} -o ${outdir}/all_samples_gene_count_id_matrix.txt \
                 -g 'gene_id' -T 12 ${bamlist} > ${logdir}/all_samples_gene_count_gene_id_matrix.log 2>&1
   featureCounts -a ${gtf} -o ${outdir}/all_samples_gene_count_name_matrix.txt \
                 -g 'gene_name' -T 12 ${bamlist} > ${logdir}/all_samples_gene_count_gene_name_matrix.log 2>&1
}
CountRepeat()
{
   indir=$1; saf=$2; outdir=$3; logdir=$4
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   bamlist=$(ls ${indir}/*.bam | tr '\n' ' '); echo ${bamlist}
   featureCounts -M -F SAF -T 12 -s 0 -a ${saf} -o ${outdir}/all_samples_repeat_count_matrix.txt \
                 ${bamlist} > ${logdir}/all_samples_repeat_count_matrix.log 2>&1
}
# TEtranscripts
TEcountPara()
{
   indir=$1; key=$2; outdir=$3; logdir=$4; ge_anno=$5; te_anno=$6; thread=$7
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   head -n 1000 ${ge_anno} > ${outdir}/gene_gtf_head1k.gtf
   ge_head1k=${outdir}/gene_gtf_head1k.gtf
   find ${indir} -name "*${key}" | xargs basename -s ".bam" \
      | parallel -P 3 --gun "TEcount -b ${indir}/{}.bam --GTF ${ge_head1k} --TE ${te_anno} --format BAM --mode multi \
                                     --sortByPos --stranded no --project ${outdir}/{}_tecount > ${logdir}/{}.log 2>&1"
}
TElocusPara()
{
   indir=$1; key=$2; outdir=$3; logdir=$4; ge_anno=$5; te_anno=$6; thread=$7
   [ ! -d "${outdir}" ] && mkdir ${outdir}; [ ! -d "${logdir}" ] && mkdir ${logdir}
   head -n 1000 ${ge_anno} > ${outdir}/gene_gtf_head1k.gtf
   ge_head1k=${outdir}/gene_gtf_head1k.gtf
   find ${indir} -name "*${key}*" | xargs basename -s ".bam" \
      | parallel -P 3 --gun "TElocus -b ${indir}/{}.bam --GTF ${ge_head1k} --TE ${te_anno} --format BAM --mode multi \
                                     --sortByPos --stranded no --project ${outdir}/{}_tecount > ${logdir}/{}.log 2>&1"
}
# Homer: make the bedgraph file to upload into UCSC
HomerBw()
{
   indir=$1; file=$2; outdir=$3; logdir=$4
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   cat ${file} | while read bam
   do
       prefix="${bam%.*}"
       makeTagDirectory ${outdir}/${prefix}_TagDir ${indir}/${bam} > ${logdir}/${prefix}_TagDir.log 2>&1
       makeUCSCfile ${outdir}/${prefix}_TagDir -o auto -fsize 5e7 -res 1 > ${logdir}/${prefix}_bedgraph.log 2>&1
       makeUCSCfile ${outdir}/${prefix}_TagDir -o auto -bigWig ${chrom} -fsize 1e20 > ${logdir}/${prefix}_bigwig.log 2>&1
   done
}
# Deeptools: create bw files
DeepCoverage(){
   indir=$1; file=$2; outdir=$3; logdir=$4; gs=$5
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   if [ ${dt} == "PE" ]; then
      cat ${file} | while read bam
      do
          prefix="${bam%.*}"
          bamCoverage --bam ${indir}/${bam} --outFileName ${outdir}/${prefix}.bw --outFileFormat bigwig \
                      --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize ${gs} \
                      --ignoreDuplicates --numberOfProcessors 16 > ${logdir}/${prefix}_coverage.log 2>&1
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


### ------------
### Run pipeline
### ------------

echo "1st step: Organize the workding directory (required)"
echo "Begin at $(date)"
cd ${wd}
mkdir logs results rawdata metadata scripts
cd logs && mkdir fastqc trimgalore multiqc star hisat2 samtools sambamba stringtie featurecounts
cd ../results && mkdir fastqc trimgalore multiqc star hisat2 samtools sambamba stringtie featurecounts && cd ../
DataLink ${rawdata} ${samples} ${wd}/rawdata
echo "Finish at $(date)"



echo "2nd step: Quality control of raw data (required)"
echo "Begin at $(date)"
ls ${wd}/rawdata/ | grep "fq.gz$" > ${wd}/metadata/raw_fq_file.txt
QcReads ${wd}/rawdata ${wd}/results/fastqc/raw ${wd}/logs/fastqc/raw
ParseQc ${wd}/results/fastqc/raw ${wd}/results/multiqc/raw ${wd}/logs/multiqc/raw
echo "Finish at $(date)"


<<'jump'
echo "3rd step: Filter and trim the reads (required)"
echo "Begin at $(date)"
quality=20; length=30
TrimReads ${wd}/rawdata ${wd}/metadata/raw_fq_file.txt ${wd}/results/trimgalore ${wd}/logs/trimgalore ${quality} ${length} ${adapter}
QcReads ${wd}/results/trimgalore ${wd}/results/fastqc/trimgalore ${wd}/logs/fastqc/trimgalore
ParseQc ${wd}/results/fastqc/trimgalore ${wd}/results/multiqc/trimgalore ${wd}/logs/multiqc/trimgalore
if [ "${dt}" == "PE" ]; then
   ls ${wd}/results/trimgalore/ | grep "val" | grep "fq.gz$" > ${wd}/metadata/trim_fq_${dt}.txt
else
   ls ${wd}/results/trimgalore/ | grep "trimmed" | grep "fq.gz$" > ${wd}/metadata/trim_fq_${dt}.txt
fi
echo "Finish at $(date)"



echo "4th step: Mapping (required)"
echo "Begin at $(date)"
# 4a. mapping with STAR (required)
for type in gene repeat
do
    # mapping
    mapping=star_${type}
    if [ "${type}" == "gene" ]; then
       StarGene ${wd}/results/trimgalore ${wd}/metadata/trim_fq_${dt}.txt ${wd}/results/star/${type} ${wd}/logs/star/${type} ${!mapping}
    elif [ "${type}" == "repeat" ]; then
       if [ "${sp}" == "mouse" ]; then
          StarRandomMus ${wd}/results/trimgalore ${wd}/metadata/trim_fq_${dt}.txt ${wd}/results/star/${type} ${wd}/logs/star/${type} ${!mapping}
       elif [ "${sp}" == "human" ]; then
          StarRandomHs ${wd}/results/trimgalore ${wd}/metadata/trim_fq_${dt}.txt ${wd}/results/star/${type} ${wd}/logs/star/${type} ${!mapping}
       else
          echo "Unrecognized species!"; exit 1
       fi
    else
       echo "Unrecognized feature types!"; exit 1
    fi
    # parse STAR results
    ParseQc ${wd}/results/star/${type} ${wd}/results/multiqc/star_${type} ${wd}/logs/multiqc/star_${type} ${type}_mapping
    if [ "${dt}" == "PE" ]; then
       ls ${wd}/results/star/${type} | grep "val" | grep "bam$" > ${wd}/metadata/star_${type}_bam_${dt}.txt
    else
       ls ${wd}/results/star/${type} | grep "trimmed" | grep "bam$" > ${wd}/metadata/star_${type}_bam_${dt}.txt
    fi
done
# 4b. mapping for repeat subfamily (optional)
<<'optional'
StarMulti ${wd}/results/trimgalore ${wd}/metadata/trim_fq_${dt}.txt ${wd}/results/star/multi ${wd}/logs/star/multi ${star_index}
ParseQc ${wd}/results/star/multi ${wd}/results/multiqc/star/multi ${wd}/logs/multiqc/star/multi multi_mapping
if [ "${dt}" == "PE" ]; then
  ls ${wd}/results/star/multi | grep "val" | grep "bam$" > ${wd}/metadata/star_multi_bam_${dt}.txt
else
  ls ${wd}/results/star/multi | grep "trimmed" | grep "bam$" > ${wd}/metadata/star_multi_bam_${dt}.txt
fi
optional
echo "Finish at $(date)"



<<'optional'
echo "5th step: Deduplication (skip this step by default)"
echo "Begin at $(date)"
BamToDedupPara ${wd}/results/star ${wd}/results/sambamba/star
if [ "${dt}" == "PE" ]; then
   ls ${wd}/results/sambamba/ | grep "val" | grep "dedup.bam$" > ${wd}/metadata/star_dedup_bam_${dt}.txt
else
   ls ${wd}/results/sambamba/ | grep "trimmed" | grep "dedup.bam$" > ${wd}/metadata/star_dedup_bam_${dt}.txt
fi
echo "Finish at $(date)"
optional



echo "8th step: Make visulization file (required)"
echo "Begin at $(date)"
# 6a. gene and repeat (required)
HomerBw ${wd}/results/star/repeat ${wd}/metadata/star_repeat_bam_${dt}.txt \
        ${wd}/results/homer/star_repeat ${wd}/logs/homer/star_repeat
# 6b. repeat subfamily (optional)
<<'optional'
HomerBw ${wd}/results/star/multi ${wd}/metadata/star_multi_bam_${dt}.txt \
        ${wd}/results/homer/star_multi ${wd}/logs/homer/star_multi ${chrom_size}
optional
echo "Finish at $(date)"



echo "9th step: Quantification (required)"
echo "Begin at $(date)"
# 7a. quantify gene and repeat (required)
CountGene ${wd}/results/star/gene ${gene_gtf} ${wd}/results/featurecounts/gene ${wd}/logs/featurecounts/gene
CountRepeat ${wd}/results/star/repeat ${repeat_saf} ${wd}/results/featurecounts/repeat ${wd}/logs/featurecounts/repeat
# 7b. repeat subfamily level (optional)
<<'optional'
TEcountPara ${wd}/results/star/multi "bam" ${wd}/results/tecount ${wd}/logs/tecount ${gene_gtf} ${repeat_gtf} 6
optional
echo "Finish at $(date)"
jump
