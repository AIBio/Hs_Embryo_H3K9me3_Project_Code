#!/bin/bash

### ---------------
### Whole procedure
### ---------------

# 1. Fastqc+MultiQC: check the quality of sequencing (required);
# 2. Trimgalore: filter and trim the reads and cut the adapter (required);
# 3. Fastqc+MultiQC: check again(required);
# 4. bwa: align the reads to genome (required);
# 5. HiCExplorer: post-process bam files (required);
# 6. FitHic: identify interactions (required);

### -------------
### Help document
### -------------

Help()
{
   echo -e ">------------------------------------------------------------------------<"
   echo "Usage: $0 -c <configure file>"
   echo "Function: This scripts was made to analyze HiC data"
   echo -e "Parameters:"
   echo -e "\t-c configure file should be stored in metadata dir:"
   echo -e "\t\t1st row: full path of working directory"
   echo -e "\t\t2nd row: full path of raw data"
   echo -e "\t\t3rd row: sample name to analyze"
   echo -e "\t\t4th row: genome index (bwa by default)"
   echo -e "\t\t5th row: chromosome size file"
   echo -e "\t\t6th row: genome sequence file (fa file)"
   echo -e "\t\t7th row: restriction enzymes, support HindIII, DpnII, BglII and MboI (format: HindIII or DpnII+BglII)"
   echo -e "\tNote: It is better to use UCSC genome file in this pipeline"
   echo -e ">------------------------------------------------------------------------<"
   exit 1
}

### ---------------
### Running command
### ---------------

# nohup bash HICv1.sh -c HICv1c.conf &

### --------------
### Parse augments
### --------------

while getopts "c:h" opt
do
   case "$opt" in
      c ) cf="$OPTARG" ;;
      h ) Help ;;
   esac
done
[ -z "${cf}" ] && Help
echo "> ---------------------------- <"
wd=`grep -v "^#" $cf | awk '{if(NR==1) print $0}'`; echo -e "[1] working directory is ${wd}"
rawdata=`grep -v "^#" $cf | awk '{if(NR==2) print $0}'`; echo -e "[2] directory of raw fastq files is ${rawdata}"
samples=`grep -v "^#" $cf | awk '{if(NR==3) print $0}'`; echo -e "[3] sample file is ${samples}"
index_dir=`grep -v "^#" $cf | awk '{if(NR==4) print $0}'`; [ ! -d ${index_dir} ] && mkdir -p ${index_dir}; echo -e "[4] genome index is ${index_dir}"
chrom_size=`grep -v "^#" $cf | awk '{if(NR==5) print $0}'`; echo -e "[5] chromosome size is ${chrom_size}"
genome_fa=`grep -v "^#" $cf | awk '{if(NR==6) print $0}'`; echo -e "[6] genome sequence file is ${genome_fa}"
enzyme=`grep -v "^#" $cf | awk '{if(NR==7) print $0}'`; echo -e "[7] restriction enzyme is ${enzyme}"
echo "> ---------------------------- <"


### ----------------
### Create functions
### ----------------
# Create soft links for raw fastq files
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
# Fastqc
QcReads(){
   indir=$1; outdir=$2; logdir=$3
   [ ! -d ${outdir} ] && mkdir -p ${outdir}; [ ! -d ${logdir} ] && mkdir -p ${logdir}
   find ${indir} -name "*.fq.gz" | xargs -P 10 -I{} fastqc {} -o ${outdir} > ${logdir}/run.log 2>&1
}
# Multiqc
ParseQc(){
   indir=$1; outdir=$2; logdir=$3
   [ ! -d ${outdir} ] && mkdir -p ${outdir}; [ ! -d ${logdir} ] && mkdir -p ${logdir}
   multiqc ${indir} -o ${outdir} > ${logdir}/run.log 2>&1
}
# Trimgalore
TrimReads(){
   indir=$1; file=$2; outdir=$3; logdir=$4
   [ ! -d ${outdir} ] && mkdir -p ${outdir}; [ ! -d ${logdir} ] && mkdir -p ${logdir}
   for i in `seq 1 2 $(cat ${file} | wc -l)`
   do
       r1=$(awk -v Num=${i} '(NR == Num){print $0}' ${file}); r1_name="${r1%%.*}"
       r2=$(awk -v Num=${i} '(NR == Num+1){print $0}' ${file}); r2_name="${r2%%.*}"
       trim_galore --paired ${indir}/${r1} ${indir}/${r2} -o ${outdir} \
                   --quality 15 --phred33 --max_n 4 --length 30 \
                   --illumina --cores 16 > ${logdir}/${r1_name}_${r2_name}.log 2>&1
   done
}
# BWA
BuildIndex(){
   seq=$1; dir=$2
   fileName="${seq##*/}"; prefix="${fileName%.*}"
   bwa index ${seq} -p ${dir}/${prefix} > ${dir}/${prefix}_run.log 2>&1
}
Mapping(){
   indir=$1; file=$2; outdir=$3; logdir=$4; index=$5
   [ ! -d ${outdir} ] && mkdir -p ${outdir}; [ ! -d ${logdir} ] && mkdir -p ${logdir}
   for i in `seq 1 2 $(cat ${file} | wc -l)`
   do
       r1=$(awk -v row=${i} '(NR == row){print $0}' ${file}); r1_name=$(basename -s ".fq.gz" ${r1})
       r2=$(awk -v row=${i} '(NR == row+1){print $0}' ${file}); r2_name=$(basename -s ".fq.gz" ${r2})
       bwa mem -t 20 -A 1 -B 4 -E 50 -L 0 ${index} ${indir}/${r1} 2>> ${logdir}/${r1_name}.log | samtools view -Shb - > ${outdir}/${r1_name}.bam
       bwa mem -t 20 -A 1 -B 4 -E 50 -L 0 ${index} ${indir}/${r2} 2>> ${logdir}/${r2_name}.log | samtools view -Shb - > ${outdir}/${r2_name}.bam
   done
}
MergeMoreBam(){
   txt=$1; outdir=$2
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}
   prefix=`echo ${txt} | xargs basename -s ".txt"`
   samtools merge -@ 16 -b ${txt} ${outdir}/${prefix}.bam
}
# HiCExplorer
BuildMX(){
   indir=$1; file=$2; outdir=$3; logdir=$4; seq=$5; enz=$6; reso=$7
   [ ! -d ${outdir} ] && mkdir -p ${outdir}/bam ${outdir}/matrix/raw ${outdir}/qc; [ ! -d ${logdir} ] && mkdir -p ${logdir}/build_matrix
   num=`echo $enz | tr "+" "\t" | awk '{print NF}'`
   if [ ${num} -eq 1 ]; then
      name=${enz}
      case "${name}" in
         "HindIII") pattern="AAGCTT" && dangling="AGCT" ;;
         "DpnII") pattern="GATC" && dangling="GATC" ;;
         "BglII") pattern="AGATCT" && dangling="GATC" ;;
         "MboI") pattern="GATC" && dangling="GATC" ;;
      esac
      cut_site=${outdir}/cut_site.bed
      [ -f ${cut_site} ] && echo "cut site file exist !" || hicFindRestSite -f ${seq} -p ${pattern} -o ${cut_site}
      for i in `seq 1 2 $(cat ${file} | wc -l)`
      do
          b1=$(awk -v row=${i} '(NR == row){print $0}' ${file}); b1_name="${b1%%.*}"
          b2=$(awk -v row=${i} '(NR == row+1){print $0}' ${file}); b2_name="${b2%%.*}"
          reso_k=`echo "${reso}/1000" | bc`; echo "Matrix resolution = ${reso}/${reso_k}kb"
          hicBuildMatrix --samFiles ${indir}/${b1} ${indir}/${b2} --threads 20 --inputBufferSize 500000 --binSize ${reso} \
                         --restrictionSequence ${pattern} --danglingSequence ${dangling} --restrictionCutFile ${cut_site} \
                         --outBam ${outdir}/bam/${b1_name}_${b2_name}.bam \
                         --outFileName ${outdir}/matrix/raw/${b1_name}_${b2_name}_${reso_k}kb.h5 \
                         --QCfolder ${outdir}/qc/${b1_name}_${b2_name} > ${logdir}/build_matrix/${b1_name}_${b2_name}_${reso_k}kb.log 2>&1
      done
   elif [ ${num} -eq 2 ]; then
      echo $enz | tr "+" "\n" | while read name
      do
          case "${name}" in
             "HindIII") pattern="AAGCTT" && dangling="AGCT" ;;
             "DpnII") pattern="GATC" && dangling="GATC" ;;
             "BglII") pattern="AGATCT" && dangling="GATC" ;;
             "MboI") pattern="GATC" && dangling="GATC" ;;
          esac
          cut_site=${outdir}/${name}_cut_site.bed
          hicFindRestSite -f ${seq} -p ${pattern} -o ${cut_site}
          echo -e "${pattern}\t${dangling}\t${cut_site}" >> ${outdir}/enzyme_anno.txt
      done
      pat1=`awk '{if(NR==1) print $1}' ${outdir}/enzyme_anno.txt`; pat2=`awk '{if(NR==2) print $1}' ${outdir}/enzyme_anno.txt`
      dang1=`awk '{if(NR==1) print $2}' ${outdir}/enzyme_anno.txt`; dang2=`awk '{if(NR==2) print $2}' ${outdir}/enzyme_anno.txt`
      cut1=`awk '{if(NR==1) print $3}' ${outdir}/enzyme_anno.txt`; cut2=`awk '{if(NR==2) print $3}' ${outdir}/enzyme_anno.txt`
      for i in `seq 1 2 $(cat ${file} | wc -l)`
      do
          b1=$(awk -v row=${i} '(NR == row){print $0}' ${file}); b1_name="${b1%%.*}"
          b2=$(awk -v row=${i} '(NR == row+1){print $0}' ${file}); b2_name="${b2%%.*}"
          hicBuildMatrix --samFiles ${indir}/${b1} ${indir}/${b2} --threads 20 --inputBufferSize 500000 --binSize 40000 \
                         --restrictionSequence ${pat1} ${pat2} --danglingSequence ${dang1} ${dang2} --restrictionCutFile ${cut1} ${cut2} \
                         --outBam ${outdir}/bam/${b1_name}_${b2_name}.bam \
                         --outFileName ${outdir}/matrix/raw/${b1_name}_${b2_name}_40kb.h5 \
                         --QCfolder ${outdir}/qc/${b1_name}_${b2_name} > ${logdir}/build_matrix/${b1_name}_${b2_name}.log 2>&1
      done
   else
      echo "The number of enzyme is higher than 2! Please modify this program and run again!"
      exit 1
   fi
}
NormalizeMX(){
   indir=$1; file=$2; outdir=$3; logdir=$4
   [ ! -d ${outdir} ] && mkdir -p ${outdir}; [ ! -d ${logdir} ] && mkdir -p ${logdir}
   in_mx_list=`awk -v path=${indir} '{print path"/"$0}' ${file} | tr "\n" " " | sed 's/ $//g'`; echo ${in_mx_list}
   out_mx_list=`awk -v path=${outdir} '{print path"/"$0}' ${file} | sed 's/\.h5/_nor.h5/g' | tr "\n" " " | sed 's/ $//g'`; echo ${out_mx_list}
   hicNormalize -m ${in_mx_list} --normalize smallest -o ${out_mx_list} > ${logdir}/run.log 2>&1
}
DiagMX(){
   indir=$1; file=$2; outdir=$3; logdir=$4
   [ ! -d ${outdir} ] && mkdir -p ${outdir}; [ ! -d ${logdir} ] && mkdir -p ${logdir}
   while read file
   do
       prefix="${file%%.*}"
       hicCorrectMatrix diagnostic_plot -m ${indir}/${file} -o ${outdir}/${prefix}_diagnostic_plot.png > ${logdir}/${prefix}_diagnostic_plot.log 2>&1
   done < ${file}
}
CorrectMx(){
   indir=$1; file=$2; outdir=$3; logdir=$4; min=$5; max=$6
   [ ! -d ${outdir} ] && mkdir -p ${outdir}; [ ! -d ${logdir} ] && mkdir -p ${logdir}
   while read file
   do
       prefix="${file%%.*}"
       # KR
       hicCorrectMatrix correct --perchr -m ${indir}/${file} -o ${outdir}/${prefix}_KR.h5 > ${logdir}/${prefix}_KR.log 2>&1
       # ICE
       hicCorrectMatrix correct --filterThreshold ${min} ${max} -m ${indir}/${file} --sequencedCountCutoff 0.2 --iterNum 500 \
                                --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
                                --correctionMethod ICE -o ${outdir}/${prefix}_ICE.h5 > ${logdir}/${prefix}_ICE.log 2>&1
   done < ${file}
}
FindTAD(){
   indir=$1; file=$2; outdir=$3; logdir=$4
   [ ! -d ${outdir} ] && mkdir -p ${outdir}; [ ! -d ${logdir} ] && mkdir -p ${logdir}
   while read matrix
   do
       prefix="${matrix%%.*}"
       hicFindTADs -m ${indir}/${matrix} --minDepth 300000 --maxDepth 3000000 --step 100000 \
                   --thresholdComparisons 0.05 --correctForMultipleTesting fdr --delta 0.01 \
                   --numberOfProcessors 20 --outPrefix ${outdir}/${prefix} > ${logdir}/${prefix}.log 2>&1
   done < ${file}
}
ABcomp(){
   indir=$1; file=$2; outdir=$3; logdir=$4
   [ ! -d ${outdir} ] && mkdir -p ${outdir}; [ ! -d ${logdir} ] && mkdir -p ${logdir}
   while read matrix
   do
       prefix="${matrix%%.*}"
       hicPCA -m ${indir}/${matrix} -o ${outdir}/${prefix}_pca1.bw ${outdir}/${prefix}_pca2.bw --format bigwig > ${logdir}/${prefix}_pca.log 2>&1
       hicTransform -m ${indir}/${matrix} -o ${outdir}/${prefix}_obs_exp.h5 --method obs_exp > ${logdir}/${prefix}_transform.log 2>&1
       hicPlotMatrix -m ${outdir}/${prefix}_pearson_all.h5 --outFileName ${outdir}/${prefix}_all_pca1.png --perChr --bigwig ${outdir}/${prefix}_all_pca1.bw
   done < ${file}
}
MergeBins(){
   indir=$1; file=$2; outdir=$3; logdir=$4; mergeN=$5
   while read matrix
   do
       prefix="${matrix%%.*}"
       hicMergeMatrixBins --matrix ${indir}/${matrix} --numBins ${mergeN} --outFileName ${outdir}/${prefix}_merge${mergeN}.h5 > ${logdir}/${prefix}_${mergeN}.log 2>&1
   done < ${file}
}
ConvFormat(){
   indir=$1; file=$2; outdir=$3; logdir=$4; inf=$5; outf=$6; res=$7
   # - indir: input directory, include the matrix file to be converted;
   # - file: file names of matrix;
   # - outdir & logdir: output matrix and log file directory;
   # - inf: input format, support hic, cool, h5, homer, HicPro;
   # - outf: output format, support cool, mcool, h5, homer, ginteractions;
   # - res: resolution of matrix, such as 10000, 20000, 40000 and so on;
   [ ! -d ${outdir} ] && mkdir -p ${outdir}; [ ! -d ${logdir} ] && mkdir -p ${logdir}
   while read mx
   do
       prefix=$(basename -s ".${inf}" ${mx})
       hicConvertFormat -m ${indir}/${mx} --inputFormat ${inf} --outputFormat ${outf} -o ${outdir}/${prefix}.${outf} --resolutions ${res}
   done < ${file}
}
TransMx(){
   # - Converts the (interaction) matrix to different types of obs/exp, pearson or covariance matrix
   indir=$1; file=$2; outdir=$3; logdir=$4; method=$5
   # - indir: input directory, include the matrix file to be converted;
   # - file: file names of matrix;
   # - outdir & logdir: output matrix and log file directory;
   # - method: transformation methods, including obs_exp, obs_exp_lieberman, obs_exp_non_zero, pearson, covariance
   [ ! -d ${outdir} ] && mkdir -p ${outdir}; [ ! -d ${logdir} ] && mkdir -p ${logdir}
   while read mx
   do
       format="${mx##*.}"
       prefix=$(basename -s ".${format}" ${mx})
       hicTransform -m ${indir}/${mx} --method ${method} -o ${outdir}/${prefix}_${method}.${format} > ${logdir}/${prefix}_${method}.log 2>&1
   done < ${file}
}
# FANC
FancPairs(){
   # - Function: Generating and filtering read Pairs
   indir=$1; file=$2; outdir=$3; logdir=$4; seq=$5; enz=$6
   # - indir: input directory, include the bam files (name-sorted);
   # - file: bam file names;
   # - outdir & logdir: output pair file and log file directory;
   # - seq: genome sequence files;
   # - enz: restriction enzyme, include HindIII, MboI;
   [ ! -d ${outdir} ] && mkdir -p ${outdir}; [ ! -d ${logdir} ] && mkdir -p ${logdir}
   for i in `seq 1 2 $(cat ${file} | wc -l)`
   do
       b1=$(awk -v row=${i} '(NR == row){print $0}' ${file}); b1n=$(basename -s ".bam" ${b1})
       b2=$(awk -v row=${i} '(NR == row+1){print $0}' ${file}); b2n=$(basename -s ".bam" ${b2})
       fanc pairs -g ${seq} -r ${enz} -m -l -p 2 -f -t 20 ${indir}/${b1} ${indir}/${b2} \
                  ${outdir}/${b1n}_${b2n}.pairs > ${logdir}/${b1n}_${b2n}.log 2>&1
       fanc pairs --re-dist-plot ${outdir}/${b1n}_${b2n}_re-dist.png ${outdir}/${b1n}_${b2n}.pairs
       fanc pairs ${outdir}/${b1n}_${b2n}.pairs --statistics-plot ${outdir}/${b1n}_${b2n}_pairs_stats.png
   done
}
FancHic(){
   # - Function: Generating, binning, and filtering Hic objects
   indir=$1; file=$2; outdir=$3; logdir=$4; reso=$5
   # - indir: input directory, include the fanc pairs files;
   # - file: pairs file names;
   # - outdir & logdir: output pair file and log file directory;
   # - reso: resolution of hic matrix, such as 10k/10000, or 1mb/1000000
   [ ! -d ${outdir} ] && mkdir -p ${outdir}; [ ! -d ${logdir} ] && mkdir -p ${logdir}
   while read pair
   do
       prefix=$(basename -s ".pairs" ${pair})
       fanc hic -n -m KR -f -t 20 -b ${reso} ${indir}/${pair} ${outdir}/${prefix}_${reso}_KR.hic > ${logdir}/${prefix}_${reso}_KR.log 2>&1
       fanc hic -n -m KR -f -t 20 -b ${reso} --low-coverage-auto ${indir}/${pair} ${outdir}/${prefix}_${reso}_KR_filter.hic > ${logdir}/${prefix}_${reso}_KR_filter.log 2>&1
       #fanc hic -n -m ICE -f -t 20 -b ${reso} ${indir}/${pair} ${outdir}/${prefix}_${reso}_ICE.hic > ${logdir}/${prefix}_${reso}_ICE.log 2>&1
   done < ${file}
}
# FitHic
FitHic(){
   indir=$1; file=$2; outdir=$3; logdir=$4; reso=$5; cs=$6
   [ ! -d ${outdir} ] && mkdir -p ${outdir}; [ ! -d ${logdir} ] && mkdir -p ${logdir}
   python /home/yhw/software/fithic-master/fithic/utils/createFitHiCFragments-fixedsize.py --chrLens ${cs} --outFile ${outdir}/fixed_fragments.txt.gz --resolution ${reso}
   cat ${file} | while read matrix
   do
       prefix="${file%%.*}"
       ConvertFormat -m ${indir}/${matrix} --inputFormat h5 --outputFormat ginteractions \
                     -o ${outdir}/${prefix}.ginteraction --resolutions ${reso}
       awk '{FS="\t"}{if($7>0) print $1,int(($3+$2)/2),$4,int(($6+$5)/2),$7}' ${outdir}/${prefix}.ginteraction.tsv | gzip - > ${outdir}/${prefix}.ginteraction.gz
       mkdir -p ${outdir}/${prefix}
       fithic -i ${outdir}/${prefix}.ginteraction.gz -f ${outdir}/fixed_fragments.txt.gz -o ${outdir}/${prefix} -r ${reso}
   done
}


### ------------
### Run pipeline
### ------------

# - 1st step: Organize the workding directory
echo "1st step: Organize the workding directory"
echo "Begin at $(date)"
cd ${wd}
mkdir logs results scripts metadata
DataLink ${rawdata} ${samples} ${wd}/rawdata
echo "Finish at $(date)"



# - 2nd step: Quality control of raw data
echo "2nd step: Quality control of raw data"
echo "Begin at $(date)"
ls ${wd}/rawdata/ | grep "fq.gz$" > ${wd}/metadata/raw_fq_file.txt
QcReads ${wd}/rawdata ${wd}/results/fastqc/raw ${wd}/logs/fastqc/raw
ParseQc ${wd}/results/fastqc/raw ${wd}/results/multiqc/raw ${wd}/logs/multiqc/raw
echo "Finish at $(date)"



# - 3rd step: Trim and filter raw reads
echo "3rd step: Filter and trim the reads"
echo "Begin at $(date)"
TrimReads ${wd}/rawdata ${wd}/metadata/raw_fq_file.txt ${wd}/results/trimgalore ${wd}/logs/trimgalore
QcReads ${wd}/results/trimgalore ${wd}/results/fastqc/trimgalore ${wd}/logs/fastqc/trimgalore
ParseQc ${wd}/results/fastqc/trimgalore ${wd}/results/multiqc/trimgalore ${wd}/logs/multiqc/trimgalore
for file in ${wd}/results/trimgalore/*fq.gz; do mv "$file" "${file/_val_2/}"; done
for file in ${wd}/results/trimgalore/*fq.gz; do mv "$file" "${file/_val_1/}"; done
ls ${wd}/results/trimgalore/ | grep "fq.gz$" > ${wd}/metadata/trim_fq_file.txt
echo "Finish at $(date)"



# - 4th step: Mapping
echo "4th step: Mapping"
# check index
[ "$(ls -A /tmp)" ] && echo "Genome index already exist" || BuildIndex ${genome_fa} ${index_dir}
tmp=`ls ${index_dir} | head -1`; prefix="${tmp%.*}"; index=${index_dir}/${prefix}
# mapping
Mapping ${wd}/results/trimgalore ${wd}/metadata/trim_fq_file.txt ${wd}/results/bwa ${wd}/logs/bwa ${index}
ls ${wd}/results/bwa/ | grep "bam$" > ${wd}/metadata/bam_file.txt
# merge replicates
for sample in 8cell_rep1 8cell_rep2 blastocyst_rep1 blastocyst_rep3 morula_rep1 morula_rep2
do
    cat metadata/bam_file.txt | grep "f1" > ${wd}/metadata/${sample}_f1_bam_merge.txt
    cat metadata/bam_file.txt | grep "r2" > ${wd}/metadata/${sample}_r2_bam_merge.txt
    awk -v path=${wd} '{print path"/results/bwa/"$0}' ${wd}/metadata/${sample}_f1_bam_merge.txt | sponge ${wd}/metadata/${sample}_f1_bam_merge.txt
    awk -v path=${wd} '{print path"/results/bwa/"$0}' ${wd}/metadata/${sample}_r2_bam_merge.txt | sponge ${wd}/metadata/${sample}_r2_bam_merge.txt
    echo "${wd}/metadata/${sample}_f1_bam_merge.txt"; cat ${wd}/metadata/${sample}_f1_bam_merge.txt
    MergeMoreBam ${wd}/metadata/${sample}_f1_bam_merge.txt ${wd}/results/samtools/merge
    echo "${wd}/metadata/${sample}_r2_bam_merge.txt"; cat ${wd}/metadata/${sample}_r2_bam_merge.txt
    MergeMoreBam ${wd}/metadata/${sample}_r2_bam_merge.txt ${wd}/results/samtools/merge
done
# sort bam files
find ${wd}/results/bwa -name "*.bam" | xargs basename -s ".bam" | xargs -P 3 -I{} samtools sort -@ 16 -n -o ${wd}/results/samtools/merge/{}_nsort.bam ${wd}/results/bwa/{}.bam
find ${wd}/results/samtools/merge -name "*.bam" | xargs basename -s ".bam" | xargs -P 3 -I{} samtools sort -@ 16 -n -o ${wd}/results/samtools/merge/{}_nsort.bam ${wd}/results/samtools/merge/{}.bam
ls ${wd}/results/samtools/merge/ | grep "nsort.bam$" > ${wd}/metadata/merged_bam_file.txt
echo "Finish at $(date)"



# - 5th step: HiCExplorer toolkit
echo "5th step: HiCExplorer toolkit"
echo "Begin at $(date)"
# 5.1 Build contact matrix (recommend a low number like 10kb)
# 10bb
BuildMX ${wd}/results/bwa ${wd}/metadata/bam_file.txt ${wd}/results/hicexplorer ${wd}/logs/hicexplorer ${genome_fa} ${enzyme} 10000
ls ${wd}/results/hicexplorer/matrix/raw/10kb | grep "CRR" | head -2 > ${wd}/metadata/raw_10kb_matrix.txt
# 40kb
BuildMX ${wd}/results/bwa ${wd}/metadata/bam_file.txt ${wd}/results/hicexplorer ${wd}/logs/hicexplorer ${genome_fa} ${enzyme} 40000
ls ${wd}/results/hicexplorer/matrix/raw/40kb | grep "40kb" > ${wd}/metadata/raw_40kb_matrix.txt

# 5.2 Normalizes matrices
# 10kb
NormalizeMX ${wd}/results/hicexplorer/matrix/raw/10kb ${wd}/metadata/raw_10kb_matrix.txt ${wd}/results/hicexplorer/matrix/nor_all_10kb ${wd}/logs/hicexplorer/matrix/nor_all_10kb
ls ${wd}/results/hicexplorer/matrix/nor_all_10kb/ | grep "10kb" > ${wd}/metadata/nor_all_10kb_matrix.txt
# 40kb
NormalizeMX ${wd}/results/hicexplorer/matrix/raw/40kb ${wd}/metadata/raw_40kb_matrix.txt ${wd}/results/hicexplorer/matrix/nor_all_40kb ${wd}/logs/hicexplorer/matrix/nor_all_40kb
ls ${wd}/results/hicexplorer/matrix/nor_all_40kb/ | grep "40kb" > ${wd}/metadata/nor_all_40kb_matrix.txt

# 5.3 Correct matrix
# 10kb
DiagMX ${wd}/results/hicexplorer/matrix/nor_all_10kb ${wd}/metadata/nor_all_10kb_matrix.txt \
       ${wd}/results/hicexplorer/matrix/correct ${wd}/logs/hicexplorer/matrix/correct
CorrectMx ${wd}/results/hicexplorer/matrix/nor_all_10kb ${wd}/metadata/nor_all_10kb_matrix.txt \
          ${wd}/results/hicexplorer/matrix/correct ${wd}/logs/hicexplorer/matrix/correct -0.5 5
ls ${wd}/results/hicexplorer/matrix/correct | grep "ICE" > ${wd}/metadata/correct_ICE_matrix.txt
# 40kb
DiagMX ${wd}/results/hicexplorer/matrix/nor_all_40kb ${wd}/metadata/nor_all_40kb_matrix.txt \
       ${wd}/results/hicexplorer/matrix/correct ${wd}/logs/hicexplorer/matrix/correct
CorrectMx ${wd}/results/hicexplorer/matrix/nor_all_40kb ${wd}/metadata/nor_all_40kb_matrix.txt \
          ${wd}/results/hicexplorer/matrix/correct ${wd}/logs/hicexplorer/matrix/correct -0.5 5
ls ${wd}/results/hicexplorer/matrix/correct | grep "ICE" > ${wd}/metadata/correct_ICE_matrix.txt

# 5.4 Convert matrix format
ConvFormat ${wd}/results/hicexplorer/matrix/correct ${wd}/metadata/correct_ICE_matrix.txt \
           ${wd}/results/hicexplorer/convformat/correct ${wd}/logs/hicexplorer/convfotmat/correct h5 cool 40000
ls ${wd}/results/hicexplorer/convformat/correct | grep "cool$" > ${wd}/metadata/correct_ICE_matrix_cool.txt

# 5.5 Transform matrix
TransMx ${wd}/results/hicexplorer/matrix/correct ${wd}/metadata/correct_ICE_matrix.txt \
        ${wd}/results/hicexplorer/transform/correct ${wd}/logs/hicexplorer/transform/correct obs_exp
ls ${wd}/results/hicexplorer/transform/correct | grep "cool$" > ${wd}/metadata/correct_ICE_matrix_obs_exp.txt

# 5.6 Identify TAD
FindTAD ${wd}/results/hicexplorer/matrix/correct ${wd}/metadata/correct_ICE_matrix.txt ${wd}/results/hicexplorer/tad ${wd}/logs/hicexplorer/tad

# 5.7 A/B compartment analysis
ABcomp ${wd}/results/hicexplorer/matrix/correct ${wd}/metadata/correct_ICE_matrix.txt ${wd}/results/hicexplorer/abs ${wd}/logs/hicexplorer/abs

# 5.8 Plotting matrix
hicPlotTADs --tracks plotTAD.ini --region chr19:56200000-58640000 -o hic_loop_zscan4_plot.png

# 5.9 Plotting aggregate
for bedpe in bedpe
do
    enhancer=./results/enhancer_${bedpe}.bed
    cut -f 1-3 ${!bedpe} > ${enhancer}
    gene=./results/gene_${bedpe}.bed
    cut -f 4-6 ${!bedpe} > ${gene}
    hicAggregateContacts --matrix ./results/hicexplorer/matrix/correct/8cell_rep1_f1_bam_merge_nsort_8cell_rep1_r2_bam_merge_nsort_40kbnor_ICE.h5 \
                         --BED ${enhancer} --BED2 ${gene} --outFileContactPairs test_tmp \
                         --outFileName 8c_rep1_aggregate_Contacts_${bedpe}_0p5to1p5.pdf --vMin 0.5 --vMax 1.5 \
                         --range 1000:20000000 --numberOfBins 100 \
                         --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
                         --avgType mean --transform obs/exp --howToCluster center --colorMap bwr 
    rm ${enhancer} ${gene}
done
echo "Finish at $(date)"



# - 6th step: FitHic: identify interactions
echo "6th step: FitHic: identify interactions"
FitHic ${wd}/results/hicexplorer/matrix/raw/40kb ${wd}/metadata/raw_40kb_matrix.txt \
       ${wd}/results/fithic ${wd}/logs/fithic 40000 ${chroms_size}
echo "Finish at $(date)"
