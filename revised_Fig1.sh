### >>> Global Setting
# - Working directory
wd=/home/cmq/bioinfo/project-cmq/embryo_93
# - H3K9me3 peaks
c4_93_pks=/home/cmq/bioinfo/project-cmq/embryo_93/results/idr/p_0.05/filtered/h4C93-Ip-1_combined_R1_val_1_h4C93-Ip-1_combined_R2_val_2_psorted_dedup_peaks_h4C93-Ip-2_combined_R1_val_1_h4C93-Ip-2_combined_R2_val_2_psorted_dedup_peaks_ucsc_bl.bed
c8_93_pks=/home/cmq/bioinfo/project-cmq/embryo_93/results/idr/p_0.05/filtered/h8C93-Ip-3_combined_R1_val_1_h8C93-Ip-3_combined_R2_val_2_psorted_dedup_peaks_h8C93-Ip-4_combined_R1_val_1_h8C93-Ip-4_combined_R2_val_2_psorted_dedup_peaks_ucsc_bl.bed
icm_93_pks=/home/cmq/bioinfo/project-cmq/embryo_93/results/idr/p_0.05/filtered/hICM93-IP-1_combined_R1_val_1_hICM93-IP-1_combined_R2_val_2_psorted_dedup_peaks_hICM93-IP-2_combined_R1_val_1_hICM93-IP-2_combined_R2_val_2_psorted_dedup_peaks_ucsc_bl.bed
morula_93_pks=/home/cmq/bioinfo/project-cmq/embryo_93/results/idr/p_0.05/filtered/hM93-Ip-1_combined_R1_val_1_hM93-Ip-1_combined_R2_val_2_psorted_dedup_peaks_hM93-Ip-2_combined_R1_val_1_hM93-Ip-2_combined_R2_val_2_psorted_dedup_peaks_ucsc_bl.bed
te_93_pks=/home/cmq/bioinfo/project-cmq/embryo_93/results/idr/p_0.05/filtered/hTE93-IP-1_combined_R1_val_1_hTE93-IP-1_combined_R2_val_2_psorted_dedup_peaks_hTE93-IP-2_combined_R1_val_1_hTE93-IP-2_combined_R2_val_2_psorted_dedup_peaks_ucsc_bl.bed
esc50_93_pks=/home/data/wanglab/HYL_low_num_ChIPseq/analysis/human/results/macs2/broad_p0.05/IP50_R1_val_1_IP50_R2_val_2_psort_dedup_6164000_bpks_ucsc_bl.bed
esc500_93_pks=/home/data/wanglab/HYL_low_num_ChIPseq/analysis/human/results/macs2/broad_p0.05/IP500_R1_val_1_IP500_R2_val_2_psort_dedup_23658746_bpks_ucsc_bl.bed
# - ATAC-seq peaks
c4_atac_pks=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/idr/dedup_bam/p_0.05/rename/c4_3pn_psorted_dedup_merged_peaks_bl_ucsc.bed
c8_atac_pks=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/idr/dedup_bam/p_0.05/rename/c8_3pn_psorted_dedup_merged_peaks_bl_ucsc.bed
icm_atac_pks=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/idr/dedup_bam/p_0.05/rename/icm_3pn_psorted_dedup_merged_peaks_bl_ucsc.bed
te_atac_pks=/home/cmq/bioinfo/project-cmq/SRP163205/results/idr/SRR7958191_peaks_SRR7958192_peaks_ucsc_bl.bed
esc_atac_pks=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/idr/dedup_bam/p_0.05/SRR5837342_SRR5837342_psorted_dedup_SRR5837343_SRR5837343_psorted_dedup_merged_peaks_SRR5837344_SRR5837344_psorted_dedup_peaks_ucsc_bl.bed
# - LiCAT-seq peaks
c4_licat_pks=/home/data/publicdata/SRP163205/analysis/results/idr/raw_0.05/rename/LiCAT_4Cell_ucsc_bl.bed
c8_licat_pks=/home/data/publicdata/SRP163205/analysis/results/idr/raw_0.05/rename/LiCAT_8Cell_ucsc_bl.bed
esc_licat_pks=/home/data/publicdata/SRP163205/analysis/results/idr/raw_0.05/rename/LiCAT_ESC_ucsc_bl.bed
icm_licat_pks=/home/data/publicdata/SRP163205/analysis/results/idr/raw_0.05/rename/LiCAT_ICM_ucsc_bl.bed
morula_licat_pks=/home/data/publicdata/SRP163205/analysis/results/idr/raw_0.05/rename/LiCAT_Morula_ucsc_bl.bed
te_licat_pks=/home/data/publicdata/SRP163205/analysis/results/idr/raw_0.05/rename/LiCAT_TE_ucsc_bl.bed
# - bed file of repeat, promoter(+-1kb around TSS) and gene body regions
re=/home/cmq/genome/ucsc/human/hg38/GRCh38_RepeatMasker_annotation_RepeatRegion_for_RegionSetEnrichment_ucsc.bed
pro1kb=/home/cmq/genome/ensembl/release97/homo_sapiens/Homo_sapiens.GRCh38.97.ens2ucsc.gene.promoter.1k.bed
pro500bp=/home/cmq/genome/ensembl/release97/homo_sapiens/Homo_sapiens.GRCh38.97.ens2ucsc.gene.promoter.500bp.bed
genebody=/home/cmq/genome/ensembl/release97/homo_sapiens/Homo_sapiens.GRCh38.97.ens2ucsc.gene.body.bed
# - Function
FunGlobalCoverageNor(){
   InDir=$1; OutDir=$2; LogDir=$3; Key=$4; Prefix=$5
   [ ! -d "${OutDir}" ] && mkdir -p ${OutDir}; [ ! -d "${LogDir}" ] && mkdir -p ${LogDir}
   BigWigList=$(ls ${InDir}/*${Key}*.bw | tr '\n' ' ')
   multiBigwigSummary bins -bs 10000 -n 0 -p 12 -b ${BigWigList} \
                           -o ${OutDir}/${Prefix}_readNor.npz \
                           --outRawCounts ${OutDir}/${Prefix}_readNor.tab > ${LogDir}/${Prefix}_nor.log 2>&1
}
FunRegionCoverageNor(){
   InDir=$1; OutDir=$2; LogDir=$3; Region=$4; Key=$5
   [ ! -d "${OutDir}" ] && mkdir -p ${OutDir}; [ ! -d "${LogDir}" ] && mkdir -p ${LogDir}
   BigWigList=$(ls ${InDir}/*${Key}*.bw | tr '\n' ' ')
   FileName="${Region##*/}"; Prefix="${FileName%.*}"
   multiBigwigSummary BED-file --BED ${Region} -p 12 -b ${BigWigList} \
                               -o ${OutDir}/${Prefix}_readNor.npz \
                               --outRawCounts ${OutDir}/${Prefix}_readNor.tab > ${LogDir}/${Prefix}_nor.log 2>&1
}



### >>> Correlation between H3K9me3 and ATAC-seq signal (CMQ)
# - Computation of genome-wide H3K9me3 signal
FunGlobalCoverageNor ${wd}/results/deeptools/bigwig/raw_bam ${wd}/results/deeptools/coverage/GobalCoverage \
                     ${wd}/logs/deeptools/coverage/GobalCoverage dedup gobal

# - Computation of H3K9me3 signal on all stages merged H3K9me3 peaks
# merge all stage peaks
cat ${wd}/results/idr/p_0.05/filtered/*ucsc_bl.bed |  sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > ${wd}/results/idr/p_0.05/filtered/covered_repeat_age/all_stage_peaks_ucsc_bl.bed
# compute signal
FunRegionCoverageNor ${wd}/results/deeptools/bigwig/raw_bam ${wd}/results/deeptools/coverage/PeaksCoverage ${wd}/logs/deeptools/coverage/PeaksCoverage \
                     ${wd}/results/idr/p_0.05/filtered/covered_repeat_age/all_stage_peaks_ucsc_bl.bed dedup

# - ATAC-seq and H3K9me3 signal on ATAC-seq and H3K9me3 merged peaks
# merge ATAC-seq and H3K9me3 peaks
cat ${wd}/results/idr/p_0.05/filtered/diff_his/k9_atac_peaks_meta.txt | while read stage pks1 pks2
do
    echo ${!pks1} ${!pks2}
    cat ${!pks1} ${!pks2} | sort -k1,1 -k2,2n | bedtools merge -i stdin \
        > ${wd}/results/idr/p_0.05/filtered/diff_his/H3K9me3_ATAC_${stage}_peaks_ucsc_bl.bed
done
# compute signal
for stage in 4Cell 8Cell ICM TE ESC50 ESC500
do 
    FunRegionCoverageNor ${wd}/results/deeptools/coverage_bigwig/atac_h3k9me3 ${wd}/results/deeptools/coverage/diff_his/h3k9me3_atac \
                         ${wd}/logs/deeptools/coverage/diff_his/h3k9me3_atac \
                         /home/cmq/bioinfo/project-cmq/embryo_93/results/idr/p_0.05/filtered/diff_his/H3K9me3_ATAC_${stage}_peaks_ucsc_bl.bed "*"
done

# - LiCAT-seq and H3K9me3 signal on LiCAT-seq and H3K9me3 merged peaks
# merge LiCAT-seq and H3K9me3 peaks
cat ${wd}/results/idr/p_0.05/filtered/diff_his/k9_licat_peaks_meta.txt | while read stage pks1 pks2
do
    echo ${!pks1} ${!pks2}
    cat ${!pks1} ${!pks2} | sort -k1,1 -k2,2n | bedtools merge -i stdin \
        > ${wd}/results/idr/p_0.05/filtered/diff_his/H3K9me3_LiCAT_${stage}_peaks_ucsc_bl.bed
done
# compute signal
for stage in 4Cell 8Cell Morula ICM TE ESC50 ESC500
do 
    FunRegionCoverageNor ${wd}/results/deeptools/coverage_bigwig/licat_h3k9me3 ${wd}/results/deeptools/coverage/diff_his/h3k9me3_licat \
                         ${wd}/logs/deeptools/coverage/diff_his/h3k9me3_licat \
                         /home/cmq/bioinfo/project-cmq/embryo_93/results/idr/p_0.05/filtered/diff_his/H3K9me3_pks/H3K9me3_${stage}_ucsc_bl.bed "*"
done



### >>> Stage specific H3K9me3 peaks annotation and heatmap visualization (YHW) 
# - stage-specific peaks were define in "revised_Fig1.R"
# - annotate specific peaks
ls *domains.bed | while read file
do
   pre=${file%%.*}
   cat ${file} | grep -v "fix" | grep -v "alt" > ${pre}_for_great.bed
done
gene=/home/cmq/genome/ensembl/release97/homo_sapiens/Homo_sapiens.GRCh38.97.ens2ucsc.gene.body.pro2kb.bed
ls *_for_great.bed | while read file
do
   name=${file%%_*}
   bedtools intersect -a ${gene} -b ${file} -f 0.1 -wa > ${name}_marked_gene_body_and_pro2kb_pcg10.txt
   bedtools intersect -a ${gene} -b ${file} -F 0.5 -wa > ${name}_marked_gene_body_and_pro2kb_gcp50.txt
done



### >>> H3K9me3 and H3K27me3 peaks enrichment analysis (CMQ)
# - H3K9me3
Rscript RegionSetEnrich.R ${wd}/results/idr/p_0.05/filtered ${re} "human" "class" 50
Rscript RegionSetEnrich.R ${wd}/results/idr/p_0.05/filtered ${pro1kb} "human" "bed" 50
Rscript RegionSetEnrich.R ${wd}/results/idr/p_0.05/filtered ${genebody} "human" "bed" 50
# - H3K27me3
Rscript RegionSetEnrich.R /home/cmq/bioinfo/project-cmq/gse124718/results/idr/multi_merged_BAMPE_p_0.05/filtered/ ${re} "human" "class" 50
Rscript RegionSetEnrich.R /home/cmq/bioinfo/project-cmq/gse124718/results/idr/multi_merged_BAMPE_p_0.05/filtered/ ${pro1kb} "human" "bed" 50
Rscript RegionSetEnrich.R /home/cmq/bioinfo/project-cmq/gse124718/results/idr/multi_merged_BAMPE_p_0.05/filtered/ ${genebody} "human" "bed" 50



### >>> H3K9me3 and H3K27me3 marked gene promoters(+-500bp around TSS) (CMQ)
# - H3K9me3 
ls ${wd}/results/idr/p_0.05/filtered/*_ucsc_bl.bed | while read pks
do 
    filename=${pks##*/}; pre=${filename%%-*}
    bedtools intersect -a ${pro500bp} -b ${pks} -f 0.25 > ${wd}/results/idr/p_0.05/filtered/peaks_covered_promoter/${pre}_covered_promoter500bp_25pect.bed
done
# - H3K27me3
ls /home/cmq/bioinfo/project-cmq/gse124718/results/idr/multi_merged_BAMPE_p_0.05/filtered/*_ucsc_bl.bed | while read pks
do 
    filename=${pks##*/}; pre=${filename%.*}
    bedtools intersect -a ${pro500bp} -b ${pks} -f 0.25 \
    > /home/cmq/bioinfo/project-cmq/gse124718/results/idr/multi_merged_BAMPE_p_0.05/filtered/peaks_covered_promoter/${pre}_covered_promoter500bp_25pect.bed
done

