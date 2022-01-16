### >>> Load functions
RegionCoverageNor(){
   indir=$1; outdir=$2; logdir=$3; key=$4; region=$5
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   bwlist=$(ls ${indir}/*${key}*.bw | tr '\n' ' ')
   fileName="${region##*/}"; prefix="${fileName%.*}"
   multiBigwigSummary BED-file --BED ${region} -p 16 -b ${bwlist} -o ${outdir}/${prefix}_readNor.npz \
                               --outRawCounts ${outdir}/${prefix}_readNor.tab > ${logdir}/${prefix}_nor.log 2>&1
}
ShiftBam(){
   indir=$1; file=$2; outdir=$3; logdir=$4
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   cat ${file} | while read bam
   do
       prefix="${bam%%.*}"
       alignmentSieve --numberOfProcessors 16 --ATACshift --bam ${indir}/${bam} -o ${outdir}/${prefix}.tmp.bam
       samtools sort -@ 16 -O bam -o ${outdir}/${prefix}_shift.bam ${outdir}/${prefix}.tmp.bam
       samtools index -@ 16 ${outdir}/${prefix}_shift.bam
       rm ${outdir}/${prefix}.tmp.bam
   done
}



### >>> Load bw files
# 4 Cell
c4_k9_fc=/home/cmq/bioinfo/project-cmq/embryo_93/results/deeptools/bigwig/merged_bam_ratio/h4C93-Ip-1_vs_h4C93-Input-1.bw
c4_atac_1=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/deeptools/bigwig/dedup_merge_bam/SRR5837324_1_val_1_SRR5837324_2_val_2_psorted_dedup_SRR5837325_1_val_1_SRR5837325_2_val_2_psorted_dedup_merged.bw
c4_atac_2=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/deeptools/bigwig/dedup_merge_bam/SRR5837326_1_val_1_SRR5837326_2_val_2_psorted_dedup_SRR5837327_1_val_1_SRR5837327_2_val_2_psorted_dedup_SRR5837328_1_val_1_SRR5837328_2_val_2_psorted_dedup_merged.bw
licat_4cell=/home/data/publicdata/SRP163205/analysis/results/deeptools/coverage/merged/SRR7958183_trimmed_psort_dedup_SRR7958184_trimmed_psort_dedup_merged.bw
# 8 Cell
c8_k9_fc=/home/cmq/bioinfo/project-cmq/embryo_93/results/deeptools/bigwig/merged_bam_ratio/h8C93-Ip-3_vs_h8C93-Input-2.bw
c8_atac_1=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/deeptools/bigwig/dedup_merge_bam/SRR5837329_1_val_1_SRR5837329_2_val_2_psorted_dedup_SRR5837330_1_val_1_SRR5837330_2_val_2_psorted_dedup_merged.bw
c8_atac_2=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/deeptools/bigwig/dedup_merge_bam/SRR5837331_1_val_1_SRR5837331_2_val_2_psorted_dedup_SRR5837332_1_val_1_SRR5837332_2_val_2_psorted_dedup_merged.bw
licat_8cell=/home/data/publicdata/SRP163205/analysis/results/deeptools/coverage/merged/SRR7958185_trimmed_psort_dedup_SRR7958186_trimmed_psort_dedup_merged.bw
# Morula
morula_k9_fc=/home/cmq/bioinfo/project-cmq/embryo_93/results/deeptools/bigwig/merged_bam_ratio/hM93-Ip-1_vs_hM93-Input.bw
morula_k273=/home/cmq/bioinfo/project-cmq/gse124718/results/deeptools/bigwig/merge_replicate_multi/morula_multi_replicate_merged.bw
licat_morula=/home/data/publicdata/SRP163205/analysis/results/deeptools/coverage/merged/SRR7958179_trimmed_psort_dedup_SRR7958180_trimmed_psort_dedup_merged.bw
# ICM
icm_k9_input=/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/hs_k9me3/results/deeptools/coverage/merge/hICM93-Input.bw
icm_k9_ip=/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/hs_k9me3/results/deeptools/coverage/merge/hICM93-IP-1.bw
icm_k9_fc=/home/cmq/bioinfo/project-cmq/embryo_93/results/deeptools/bigwig/merged_bam_ratio/hICM93-IP-1_vs_hICM93-Input.bw
icm_k273=/home/cmq/bioinfo/project-cmq/gse124718/results/deeptools/bigwig/merge_replicate_multi/icm_multi_replicate_merged.bw
icm_k43=/home/cmq/bioinfo/project-cmq/gse124718_human_43/results/deeptools/bigwig/merged_bam/SRR8402635_trimmed_psorted_dedup_SRR8402636_trimmed_psorted_dedup_merged.bw
icm_27ac=/home/yhw/bioinfo/project-temp/wjc/gse124718/results/deeptools/coverage/SRR9131737_1_val_1_SRR9131737_2_val_2_psort_dedup.bw
icm_atac_1=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/deeptools/bigwig/dedup_merge_bam/SRR5837337_1_val_1_SRR5837337_2_val_2_psorted_dedup_SRR5837338_1_val_1_SRR5837338_2_val_2_psorted_dedup_merged.bw
icm_atac_2=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/deeptools/bigwig/dedup_merge_bam/SRR5837339_1_val_1_SRR5837339_2_val_2_psorted_dedup_SRR5837340_1_val_1_SRR5837340_2_val_2_psorted_dedup_SRR5837341_1_val_1_SRR5837341_2_val_2_psorted_dedup_merged.bw
icm_atac_3pn1=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/deeptools/bigwig/dedup_merge_bam/rename/icm_3pn_1_dedup_merged.bw
icm_atac_3pn2=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/deeptools/bigwig/dedup_merge_bam/rename/icm_3pn_2_dedup_merged.bw
licat_icm=/home/data/publicdata/SRP163205/analysis/results/deeptools/coverage/merged/SRR7958181_trimmed_psort_dedup_SRR7958182_trimmed_psort_dedup_merged.bw
# TE
te_k9_input=/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/hs_k9me3/results/deeptools/coverage/merge/hTE93-Input.bw
te_k9_ip=/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/hs_k9me3/results/deeptools/coverage/merge/hTE93-IP-1.bw
te_k9_fc=/home/cmq/bioinfo/project-cmq/embryo_93/results/deeptools/bigwig/merged_bam_ratio/hTE93-IP-1_vs_hTE93-Input.bw
te_k273=/home/cmq/bioinfo/project-cmq/gse124718/results/deeptools/bigwig/merge_replicate_multi/te_multi_replicate_merged.bw
licat_te=/home/data/publicdata/SRP163205/analysis/results/deeptools/coverage/merged/SRR7958191_trimmed_psort_dedup_SRR7958192_trimmed_psort_dedup_merged.bw
# Naive & Primed
naive_93_fc1=/home/data/publicdata/GSE84382/analysis/results/deeptools/coverage/foldchange/SRR3904795_trimmed_psort_dedup_vs_SRR3904791_trimmed_psort_dedup.bw
naive_93_fc2=/home/data/publicdata/GSE84382/analysis/results/deeptools/coverage/foldchange/SRR3904797_trimmed_psort_dedup_vs_SRR3904791_trimmed_psort_dedup.bw
naive_273_fc=/home/data/publicdata/GSE59434/analysis/results/deeptools/bamcompare/dedup_bam/WIBR2_6iLA_H3K27me3_trimmed_psort_dedup_vs_WIBR2_6iLA_Input_trimmed_psort_dedup.bw
naive_43_fc=/home/data/publicdata/GSE59434/analysis/results/deeptools/bamcompare/dedup_bam/WIBR2_6iLA_H3K4me3_trimmed_psort_dedup_vs_WIBR2_6iLA_Input_trimmed_psort_dedup.bw
naive_27ac_fc=/home/data/publicdata/GSE69646/analysis/results/deeptools/bamcompare/Naive_H3K27ac_rep1_trimmed_psort_dedup_vs_Naive_Input_trimmed_psort_dedup.bw
primed_93_fc1=/home/data/publicdata/GSE84382/analysis/results/deeptools/coverage/foldchange/SRR3904792_trimmed_psort_dedup_vs_SRR3904789_trimmed_psort_dedup.bw
primed_93_fc2=/home/data/publicdata/GSE84382/analysis/results/deeptools/coverage/foldchange/SRR3904794_trimmed_psort_dedup_vs_SRR3904789_trimmed_psort_dedup.bw
primed_273_fc=/home/data/publicdata/GSE59434/analysis/results/deeptools/bamcompare/dedup_bam/WIBR2_hESM_H3K27me3_trimmed_psort_dedup_vs_WIBR2_hESM_Input_trimmed_psort_dedup.bw
primed_43_fc=/home/data/publicdata/GSE59434/analysis/results/deeptools/bamcompare/dedup_bam/WIBR2_hESM_H3K4me3_trimmed_psort_dedup_vs_WIBR2_hESM_Input_trimmed_psort_dedup.bw
primed_27ac_fc=/home/data/publicdata/GSE69646/analysis/results/deeptools/bamcompare/Primed_H3K27ac_rep1_trimmed_psort_dedup_vs_Primed_Input_trimmed_psort_dedup.bw
# Normalized H3K27me3
icm_k273=/home/cmq/bioinfo/project-cmq/gse124718/results/deeptools/bigwig/merge_replicate_multi_scale/icm_multi_replicate_merged.bw
morula_k273=/home/cmq/bioinfo/project-cmq/gse124718/results/deeptools/bigwig/merge_replicate_multi_scale/morula_multi_replicate_merged.bw
te_k273=/home/cmq/bioinfo/project-cmq/gse124718/results/deeptools/bigwig/merge_replicate_multi_scale/te_multi_replicate_merged.bw
# New TE K9 LFC
te_k93_1=./results/deeptools/bamcompare/new/hTE93-IP-1_combined_R1_val_1_hTE93-IP-1_combined_R2_val_2_psorted_dedup_hTE93-IP-2_combined_R1_val_1_hTE93-IP-2_combined_R2_val_2_psorted_dedup_merged_vs_h4C93-Input-1_combined_R1_val_1_h4C93-Input-1_combined_R2_val_2_psorted_dedup.bw
te_k93_2=./results/deeptools/bamcompare/new/hTE93-IP-1_combined_R1_val_1_hTE93-IP-1_combined_R2_val_2_psorted_dedup_hTE93-IP-2_combined_R1_val_1_hTE93-IP-2_combined_R2_val_2_psorted_dedup_merged_vs_h8C93-Input-2_combined_R1_val_1_h8C93-Input-2_combined_R2_val_2_psorted_dedup.bw
te_k93_3=./results/deeptools/bamcompare/new/hTE93-IP-1_combined_R1_val_1_hTE93-IP-1_combined_R2_val_2_psorted_dedup_hTE93-IP-2_combined_R1_val_1_hTE93-IP-2_combined_R2_val_2_psorted_dedup_merged_vs_hICM93-Input_combined_R1_val_1_hICM93-Input_combined_R2_val_2_psorted_dedup.bw
te_k93_4=./results/deeptools/bamcompare/new/hTE93-IP-1_combined_R1_val_1_hTE93-IP-1_combined_R2_val_2_psorted_dedup_hTE93-IP-2_combined_R1_val_1_hTE93-IP-2_combined_R2_val_2_psorted_dedup_merged_vs_hM93-Input_combined_R1_val_1_hM93-Input_combined_R2_val_2_psorted_dedup.bw
te_k93_5=./results/deeptools/bamcompare/new/hTE93-IP-1_combined_R1_val_1_hTE93-IP-1_combined_R2_val_2_psorted_dedup_hTE93-IP-2_combined_R1_val_1_hTE93-IP-2_combined_R2_val_2_psorted_dedup_merged_vs_hTE93-Input_combined_R1_val_1_hTE93-Input_combined_R2_val_2_psorted_dedup.bw



### >>> Load peaks
# Morula
morula_k9_pks=/home/yhw/bioinfo/project-mine/Embryo.93/R/RawData/peaks/hs_k9/hM93-Ip-1_combined_R1_val_1_hM93-Ip-1_combined_R2_val_2_psorted_dedup_peaks_hM93-Ip-2_combined_R1_val_1_hM93-Ip-2_combined_R2_val_2_psorted_dedup_peaks_ucsc_bl.bed
morula_k27_pks=/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/gse123023/results/idr/broad_p0.05/filtered/SRR8249175_1_val_1_SRR8249175_2_val_2_psort_dedup_peaks_SRR8249176_1_val_1_SRR8249176_2_val_2_psort_dedup_peaks_ucsc_bl.bed
# ICM
icm_k9_pks=/home/yhw/bioinfo/project-mine/Embryo.93/R/RawData/peaks/hs_k9/hICM93-IP-1_combined_R1_val_1_hICM93-IP-1_combined_R2_val_2_psorted_dedup_peaks_hICM93-IP-2_combined_R1_val_1_hICM93-IP-2_combined_R2_val_2_psorted_dedup_peaks_ucsc_bl.bed
icm_k27_pks=/home/cmq/bioinfo/project-cmq/gse124718/results/idr/multi_merged_BAMPE_p_0.01/SRR8402638_trimmed_psorted_dedup_peaks_SRR8402639_trimmed_psorted_dedup_peaks_ucsc_bl.bed
icm_atac_1_pks=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/idr/dedup_bam/p_0.05/SRR5837337_SRR5837337_psorted_dedup_SRR5837338_SRR5837338_psorted_dedup_merged_peaks_SRR5837339_SRR5837339_psorted_dedup_SRR5837340_SRR5837340_psorted_dedup_SRR5837341_SRR5837341_psorted_dedup_merged_peaks_ucsc.bed
icm_k27ac_pks=/home/yhw/bioinfo/project-temp/wjc/gse124718/results/idr/narrow_p0.05/SRR9131737_1_val_1_SRR9131737_2_val_2_psort_dedup_peaks_SRR9131738_1_val_1_SRR9131738_2_val_2_psort_dedup_peaks_bl.bed
licat_icm_pks=/home/data/publicdata/SRP163205/analysis/results/idr/raw_0.05/SRR7958181_trimmed_psort_dedup_peaks_SRR7958182_trimmed_psort_dedup_peaks_ucsc_bl.bed
# TE
te_k9_pks=/home/yhw/bioinfo/project-mine/Embryo.93/R/RawData/peaks/hs_k9/hTE93-IP-1_combined_R1_val_1_hTE93-IP-1_combined_R2_val_2_psorted_dedup_peaks_hTE93-IP-2_combined_R1_val_1_hTE93-IP-2_combined_R2_val_2_psorted_dedup_peaks_ucsc_bl.bed
te_k27_pks=/home/cmq/bioinfo/project-cmq/gse124718/results/idr/multi_merged_BAMPE_p_0.05/filtered/SRR9131731_SRR9131731_psorted_dedup_peaks_SRR9131732_SRR9131732_psorted_dedup_peaks_ucsc_bl.bed
licat_te_pks=/home/data/publicdata/SRP163205/analysis/results/idr/raw_0.05/SRR7958191_trimmed_psort_dedup_peaks_SRR7958192_trimmed_psort_dedup_peaks_ucsc_bl.bed
# promoter
hs_pro3k=/home/yhw/genome/ensembl/release97/homo_sapiens/dna_anno/regions/Homo_sapiens.GRCh38.97_pro3k_ucsc.bed
hs_pro1k=/home/yhw/genome/ensembl/release97/homo_sapiens/dna_anno/regions/Homo_sapiens.GRCh38.97_pro1k_ucsc.bed
# morula single cell
morula_outer_open=/home/yhw/bioinfo/project-mine/Embryo.93/scATAC/Morula_outer_cell_specific_open_chromatin_fc1.5_hg38.bed
morula_inner_open=/home/yhw/bioinfo/project-mine/Embryo.93/scATAC/Morula_inner_cell_specific_open_chromatin_fc1.5_hg38.bed
# ICM & TE specific gained peaks
icm_spe_gain_93_pks=/home/cmq/bioinfo/project-cmq/embryo_93/results/idr/p_0.05/filtered/morula_icm_te/morula_to_icm_gained_specific_peaks_compared_to_te_ucsc_bl.bed
te_spe_gain_93_pks=/home/cmq/bioinfo/project-cmq/embryo_93/results/idr/p_0.05/filtered/morula_icm_te/morula_to_te_gained_specific_peaks_compared_to_icm_ucsc_bl.bed
# filter ICM & TE specific gained H3K9me3 peaks by ATAC-seq peaks
bedtools intersect -a ${icm_spe_gain_93_pks} -b ${icm_atac_1_pks} -v > ./results/bedtools/icm_vs_te/morula_to_icm_gained_specific_peaks_compared_to_te_ucsc_bl_filter_by_atac.bed
bedtools intersect -a ${te_spe_gain_93_pks} -b ${licat_te_pks} -v > ./results/bedtools/icm_vs_te/morula_to_te_gained_specific_peaks_compared_to_icm_ucsc_bl_filter_by_licat.bed
icm_spe_gain_93_pks=/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/hs_k9me3/results/bedtools/icm_vs_te/morula_to_icm_gained_specific_peaks_compared_to_te_ucsc_bl_filter_by_atac.bed
te_spe_gain_93_pks=/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/hs_k9me3/results/bedtools/icm_vs_te/morula_to_te_gained_specific_peaks_compared_to_icm_ucsc_bl_filter_by_licat.bed
# filter peaks with morula ATAC peaks
morula_atac_pks=/home/yhw/bioinfo/project-mine/Embryo.93/scATAC/Morula_all_open_chromatin_peaks.bed
bedtools intersect -a ${icm_spe_gain_93_pks} -b ${morula_atac_pks} -F 0.5 -wa -wb > ./results/bedtools/icm_vs_te/morula_to_icm_gained_specific_peaks_compared_to_te_ucsc_bl_filter_by_atac_in_morula.bed
bedtools intersect -a ${te_spe_gain_93_pks} -b ${morula_atac_pks} -F 0.5 -wa -wb > ./results/bedtools/icm_vs_te/morula_to_te_gained_specific_peaks_compared_to_icm_ucsc_bl_filter_by_licat_in_morula.bed



### >>> ICM and TE specific gained K9 peaks in morula (YHW)
# open in opposite lineage
icm_spe_gain_93_pks_open_in_TE=${wd}/results/bedtools/icm_vs_te/morula_to_icm_gained_specific_peaks_compared_to_te_ucsc_bl_filter_by_atac_open_in_TE.bed
te_spe_gain_93_pks_open_in_ICM=${wd}/results/bedtools/icm_vs_te/morula_to_te_gained_specific_peaks_compared_to_icm_ucsc_bl_filter_by_licat_open_in_ICM.bed
bedtools intersect -a ${icm_spe_gain_93_pks} -b ${licat_te_pks} | awk '{if(($3-$2)>=200) print $0}' > ${icm_spe_gain_93_pks_open_in_TE}
bedtools intersect -a ${te_spe_gain_93_pks} -b ${icm_atac_1_pks} | awk '{if(($3-$2)>=200) print $0}' > ${te_spe_gain_93_pks_open_in_ICM}
cat ${icm_spe_gain_93_pks_open_in_TE} ${te_spe_gain_93_pks_open_in_ICM} | cut -f 1-3 | sort -k1,1 -k2,2n > ${wd}/results/bedtools/icm_vs_te/morula_to_icm_te_gained_specific_peaks_filter_by_atac_open.bed
# close in opposite lineage
icm_spe_gain_93_pks_close_in_TE=${wd}/results/bedtools/icm_vs_te/morula_to_icm_gained_specific_peaks_compared_to_te_ucsc_bl_filter_by_atac_close_in_TE.bed
te_spe_gain_93_pks_close_in_ICM=${wd}/results/bedtools/icm_vs_te/morula_to_te_gained_specific_peaks_compared_to_icm_ucsc_bl_filter_by_licat_close_in_ICM.bed
bedtools intersect -a ${icm_spe_gain_93_pks} -b ${licat_te_pks} -v -wa > ${icm_spe_gain_93_pks_close_in_TE}
bedtools intersect -a ${te_spe_gain_93_pks} -b ${icm_atac_1_pks} -v -wa > ${te_spe_gain_93_pks_close_in_ICM}



### >>> Inner and outer specific open chromatin in ICM and TE (YHW)
# H3K9me3 binds to the open chromatin
morula_outer_open_icm_k9=/home/yhw/bioinfo/project-mine/Embryo.93/scATAC/Morula_outer_cell_specific_open_chromatin_fc1.5_hg38_intersect_with_icm_k9.bed
morula_inner_open_te_k9=/home/yhw/bioinfo/project-mine/Embryo.93/scATAC/Morula_inner_cell_specific_open_chromatin_fc1.5_hg38_intersect_with_te_k9.bed
bedtools intersect -a ${morula_outer_open} -b ${icm_k9_pks} | sort -k1,1 -k2,2n > ${morula_outer_open_icm_k9}
bedtools intersect -a ${morula_inner_open} -b ${te_k9_pks} | sort -k1,1 -k2,2n > ${morula_inner_open_te_k9}
# H3K27me3 binds to the open chromatin
morula_outer_open_icm_k27=/home/yhw/bioinfo/project-mine/Embryo.93/scATAC/Morula_outer_cell_specific_open_chromatin_fc1.5_hg38_intersect_with_icm_k27.bed
morula_inner_open_te_k27=/home/yhw/bioinfo/project-mine/Embryo.93/scATAC/Morula_inner_cell_specific_open_chromatin_fc1.5_hg38_intersect_with_te_k27.bed
bedtools intersect -a ${morula_outer_open} -b ${icm_k27_pks} | sort -k1,1 -k2,2n > ${morula_outer_open_icm_k27}
bedtools intersect -a ${morula_inner_open} -b ${te_k27_pks} | sort -k1,1 -k2,2n > ${morula_inner_open_te_k27}



### >>> Plotting (YHW)
# morula to blastocyst
for pks in morula_outer_open_icm_k27 morula_inner_open_te_k27 morula_outer_open_icm_k9 morula_inner_open_te_k9 #morula_outer_open morula_inner_open
do
    mkdir -p ./results/deeptools/matrix/reference/icm_vs_te; mkdir -p ./results/deeptools/heatmap/reference/icm_vs_te
    mkdir -p ./results/deeptools/profile/reference/icm_vs_te
    # H3K9me3
    computeMatrix reference-point -S ${morula_k273} ${icm_k273} ${te_k273} \
                                  -R ${!pks} \
                                  --averageTypeBins mean --referencePoint center \
                                  --binSize 50 --beforeRegionStartLength 3000 --afterRegionStartLength 3000 \
                                  --numberOfProcessors 40 --missingDataAsZero -o ./results/deeptools/matrix/reference/icm_vs_te/${pks}_K273_RPKM.matrix.gz
    plotProfile -m ./results/deeptools/matrix/reference/icm_vs_te/${pks}_K273_RPKM.matrix.gz \
                -out ./results/deeptools/profile/reference/icm_vs_te/${pks}_K273_RPKM.matrix.pdf \
                --plotType=lines --perGroup --plotHeight 6 --plotWidth 6 --samplesLabel "Morula 273" "ICM 273" "TE 273" --yMin 0 --yMax 25
    # H3K27me3
    computeMatrix reference-point -S ${morula_k9_fc} ${icm_k9_fc} ${te_k9_fc} \
                                  -R ${!pks} \
                                  --averageTypeBins mean --referencePoint center \
                                  --binSize 50 --beforeRegionStartLength 3000 --afterRegionStartLength 3000 \
                                  --numberOfProcessors 40 --missingDataAsZero -o ./results/deeptools/matrix/reference/icm_vs_te/${pks}_K9_LFC.matrix.gz
    plotProfile -m ./results/deeptools/matrix/reference/icm_vs_te/${pks}_K9_LFC.matrix.gz \
                -out ./results/deeptools/profile/reference/icm_vs_te/${pks}_K9_LFC.matrix.pdf \
                --plotType=lines --perGroup --plotHeight 6 --plotWidth 6 --samplesLabel "Morula 93" "ICM 93" "TE 93" --yMin 0 --yMax 3
done
# ICM and TE specific peaks
mkdir -p ./results/deeptools/matrix/reference/icm_vs_te; mkdir -p ./results/deeptools/heatmap/reference/icm_vs_te
computeMatrix reference-point -S ${c4_k9_fc} ${c8_k9_fc} ${morula_k9_fc} ${icm_k9_fc} ${te_k9_fc} \
                              -R ${icm_spe_gain_93_pks} ${te_spe_gain_93_pks} --averageTypeBins mean --referencePoint center \
                              --binSize 50 --beforeRegionStartLength 3000 --afterRegionStartLength 3000 \
                              --numberOfProcessors 40 --missingDataAsZero -o ./results/deeptools/matrix/reference/icm_vs_te/icm_te_spe_gain_93_pks_K93_LFC.matrix.gz
plotHeatmap -m ./results/deeptools/matrix/reference/icm_vs_te/icm_te_spe_gain_93_pks_K93_LFC.matrix.gz \
            -o ./results/deeptools/heatmap/reference/icm_vs_te/icm_te_spe_gain_93_pks_K93_LFC.matrix.pdf \
            --averageTypeSummaryPlot mean --colorList "cyan,white,red" --heatmapHeight 6 --heatmapWidth 2 --whatToShow "plot, heatmap and colorbar" \
            --plotType lines --zMin -1 --zMax 3 --legendLocation upper-left --interpolationMethod auto \
            --samplesLabel "4 Cell" "8 Cell" "Morula" "ICM" "TE"
# abnormal embryo
abnormal_down=./metadata/Venn_diagram_DERs_between_abembryo_8cell_746.bed
mkdir -p ./results/deeptools/matrix/reference/abnormal_embryo; mkdir -p ./results/deeptools/heatmap/reference/abnormal_embryo
computeMatrix reference-point -S ${c4_k9_fc} ${c8_k9_fc} ${morula_k9_fc} ${icm_k9_fc} ${te_k9_fc} \
                              -R ${abnormal_down} --averageTypeBins mean --referencePoint center \
                              --binSize 50 --beforeRegionStartLength 3000 --afterRegionStartLength 3000 \
                              --numberOfProcessors 40 --missingDataAsZero -o ./results/deeptools/matrix/reference/icm_vs_te/icm_te_spe_gain_93_pks_K93_LFC.matrix.gz
plotHeatmap -m ./results/deeptools/matrix/reference/icm_vs_te/icm_te_spe_gain_93_pks_K93_LFC.matrix.gz \
            -o ./results/deeptools/heatmap/reference/icm_vs_te/icm_te_spe_gain_93_pks_K93_LFC.matrix.pdf \
            --averageTypeSummaryPlot mean --colorList "cyan,white,red" --heatmapHeight 6 --heatmapWidth 2 --whatToShow "plot, heatmap and colorbar" \
            --plotType lines --zMin -1 --zMax 3 --legendLocation upper-left --interpolationMethod auto \
            --samplesLabel "4 Cell" "8 Cell" "Morula" "ICM" "TE"



### >>> Transposons marked by H3K9me3 in TE (YHW)
repeat_by_te_k9=/home/cmq/bioinfo/project-cmq/embryo_93/results/idr/p_0.05/filtered/morula_icm_te/enrich/morula_to_te_gained_specific_peaks_compared_to_icm_ucsc_bl_filter_by_licat_marked_repeats_subfamily_num20_rat10_region.bed
# classification
repeat_by_te_k9_icm_atac=./results/bedtools/icm_vs_te/morula_to_te_gained_specific_peaks_marked_repeats_in_ICM_ATAC.bed
repeat_by_te_k9_icm_27ac=./results/bedtools/icm_vs_te/morula_to_te_gained_specific_peaks_marked_repeats_in_ICM_H3K27ac.bed
repeat_by_te_k9_icm_27m3=./results/bedtools/icm_vs_te/morula_to_te_gained_specific_peaks_marked_repeats_in_ICM_H3K27me3.bed
bedtools intersect -a ${repeat_by_te_k9} -b ${icm_atac_1_pks} | sort -k1,1 -k2,2n > ${repeat_by_te_k9_icm_atac}
bedtools intersect -a ${repeat_by_te_k9} -b ${icm_k27ac_pks} | sort -k1,1 -k2,2n > ${repeat_by_te_k9_icm_27ac}
bedtools intersect -a ${repeat_by_te_k9} -b ${icm_k27_pks} | sort -k1,1 -k2,2n > ${repeat_by_te_k9_icm_27m3}
cat ${repeat_by_te_k9_icm_atac} ${repeat_by_te_k9_icm_27ac} | sort -k1,1 -k2,2n > results/bedtools/icm_vs_te/morula_to_te_gained_specific_peaks_marked_repeats_in_ICM_ATAC_H3K27ac_filtered.bed
bedtools intersect -v -a results/bedtools/icm_vs_te/morula_to_te_gained_specific_peaks_marked_repeats_in_ICM_ATAC_H3K27ac.bed -b ${licat_te_pks} > results/bedtools/icm_vs_te/morula_to_te_gained_specific_peaks_marked_repeats_in_ICM_ATAC_H3K27ac_filtered2.bed
# compute RPKM values (regarded as the score column)
RegionCoverageNor ./results/deeptools/coverage/icm_vs_te ./results/deeptools/multibw/icm_vs_te ./logs/deeptools/multibw/icm_vs_te "*" ${repeat_by_te_k9_icm_atac_1}
cat results/deeptools/multibw/icm_vs_te/morula_to_te_gained_specific_peaks_marked_repeats_in_ICM_ATAC_readNor.tab \
   | grep -v "_" | sort -k1,1 -k2,2n \
   | awk '{print $0"\t""*"}' > results/deeptools/multibw/icm_vs_te/morula_to_te_gained_specific_peaks_marked_repeats_in_ICM_ATAC_readNor.bed

# plotting
mkdir -p ./results/deeptools/matrix/reference/icm_vs_te ./results/deeptools/heatmap/reference/icm_vs_te
computeMatrix reference-point -S ${icm_27ac} ${icm_atac_1} ${icm_atac_2} ${icm_k273} \
                              -R ${repeat_by_te_k9_icm_atac} ${repeat_by_te_k9_icm_27ac} ${repeat_by_te_k9_icm_27m3} \
                              --averageTypeBins mean --referencePoint center \
                              --binSize 50 --beforeRegionStartLength 3000 --afterRegionStartLength 3000 \
                              --numberOfProcessors 40 --missingDataAsZero \
                              -o ./results/deeptools/matrix/reference/icm_vs_te/repeat_by_te_k9_icm_atac_27ac_27me3_RPKM.matrix.gz
for min in 0
do
    for max in `seq 30 5 40`
    do
        for color in Blues Greens Oranges
        do
            plotHeatmap -m ./results/deeptools/matrix/reference/icm_vs_te/repeat_by_te_k9_icm_atac_27ac_27me3_RPKM.matrix.gz \
                        -o ./results/deeptools/heatmap/reference/icm_vs_te/repeat_by_te_k9_icm_atac_27ac_27me3_RPKM_${min}-${max}_${color}.matrix.pdf \
                        --averageTypeSummaryPlot mean --colorMap ${color} --heatmapHeight 6 --heatmapWidth 2 --whatToShow "plot, heatmap and colorbar" \
                        --plotType lines --zMin ${min} --zMax ${max} --legendLocation upper-left --interpolationMethod auto \
                        --samplesLabel "H3K27ac" "ATAC.1" "ATAC.2" "H3K27me3"
        done
    done
done



### >>> Footprint analysis (YHW)
# ATAC-seq processing from TE --- Liu, L. et al.  Nat Commun (2019)
te_atac_bam_dir=/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/hs_k9me3/results/samtools/te_atac
ls ${te_atac_bam_dir} | grep "bam$" > ${te_atac_bam_dir}/bam_files.txt
ShiftBam ${te_atac_bam_dir} ${te_atac_bam_dir}/bam_files.txt ${te_atac_bam_dir} ${te_atac_bam_dir}
# footprint
genome_fa=/home/yhw/genome/ucsc/homo_sapiens/hg38/hg38.fa
blacklist=/home/yhw/software/blender/hg38.blacklist.bed
icm_bam=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/samtools/merge/dedup_bam/SRR5837337_1_val_1_SRR5837337_2_val_2_psorted_dedup_SRR5837338_1_val_1_SRR5837338_2_val_2_psorted_dedup_merged.bam
te_bam=/home/data/publicdata/SRP163205/analysis/results/samtools/merge/SRR7958191_trimmed_psort_dedup_SRR7958192_trimmed_psort_dedup_merged.bam
primed_bam1=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/samtools/merge/dedup_bam/SRR5837344_1_val_1_SRR5837344_2_val_2_psorted_dedup.bam
primed_bam2=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/samtools/merge/dedup_bam/SRR5837342_1_val_1_SRR5837342_2_val_2_psorted_dedup_SRR5837343_1_val_1_SRR5837343_2_val_2_psorted_dedup_merged.bam
naive_bam1=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/samtools/merge/dedup_bam/SRR5878328_1_val_1_SRR5878328_2_val_2_psorted_dedup_SRR5878329_1_val_1_SRR5878329_2_val_2_psorted_dedup_merged.bam
naive_bam2=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/samtools/merge/dedup_bam/SRR6961739_1_val_1_SRR6961739_2_val_2_psorted_dedup_SRR5878330_1_val_1_SRR5878330_2_val_2_psorted_dedup_merged.bam
binding_pks=./results/bedtools/icm_vs_te/morula_to_te_gained_specific_peaks_marked_repeats_in_ICM_ATAC_H3K27ac_filtered.bed
jaspar=/home/yhw/document/jaspar/v2020/jaspar2020.jaspar
mkdir -p ./results/tobias/ICM_ATAC_H3K27ac_filter3
binding_pks=${wd}/results/bedtools/icm_vs_te/morula_to_icm_te_gained_specific_peaks_filter_by_atac_open.bed
jaspar=${wd}/metadata/motif.jaspar
# ATACorrect
TOBIAS ATACorrect --bam ${icm_bam} --genome ${genome_fa} --peaks ${binding_pks} --blacklist ${blacklist} --outdir ./results/tobias/ICM_ATAC_H3K27ac_filter3 --cores 20 --prefix ICM
TOBIAS ATACorrect --bam ${te_bam} --genome ${genome_fa} --peaks ${binding_pks} --blacklist ${blacklist} --outdir ./results/tobias/ICM_ATAC_H3K27ac_filter3 --cores 20 --prefix TE
TOBIAS ATACorrect --bam ${naive_bam1} --genome ${genome_fa} --peaks ${binding_pks} --blacklist ${blacklist} --outdir ./results/tobias/ICM_ATAC_H3K27ac_filter3 --cores 20 --prefix Naive1
TOBIAS ATACorrect --bam ${naive_bam2} --genome ${genome_fa} --peaks ${binding_pks} --blacklist ${blacklist} --outdir ./results/tobias/ICM_ATAC_H3K27ac_filter3 --cores 20 --prefix Naive2
TOBIAS ATACorrect --bam ${primed_bam1} --genome ${genome_fa} --peaks ${binding_pks} --blacklist ${blacklist} --outdir ./results/tobias/ICM_ATAC_H3K27ac_filter3 --cores 20 --prefix Primed1
TOBIAS ATACorrect --bam ${primed_bam2} --genome ${genome_fa} --peaks ${binding_pks} --blacklist ${blacklist} --outdir ./results/tobias/ICM_ATAC_H3K27ac_filter3 --cores 20 --prefix Primed2
icm_cbw=${wd}/results/tobias/ICM_ATAC_H3K27ac_filter3/ICM_corrected.bw
te_cbw=${wd}/results/tobias/ICM_ATAC_H3K27ac_filter3/TE_corrected.bw
naive_cbw1=${wd}/results/tobias/ICM_ATAC_H3K27ac_filter3/Naive1_corrected.bw
naive_cbw2=${wd}/results/tobias/ICM_ATAC_H3K27ac_filter3/Naive2_corrected.bw
primed_cbw1=${wd}/results/tobias/ICM_ATAC_H3K27ac_filter3/Primed1_corrected.bw
primed_cbw2=${wd}/results/tobias/ICM_ATAC_H3K27ac_filter3/Primed2_corrected.bw
# ScoreBigwig
TOBIAS ScoreBigwig --signal ${icm_cbw} --regions ${binding_pks} --output ${wd}/results/tobias/ICM_ATAC_H3K27ac_filter3/ICM_footprints.bw --cores 20
TOBIAS ScoreBigwig --signal ${te_cbw} --regions ${binding_pks} --output ${wd}/results/tobias/ICM_ATAC_H3K27ac_filter3/TE_footprints.bw --cores 20
icm_fbw=${wd}/results/tobias/ICM_ATAC_H3K27ac_filter3/ICM_footprints.bw
te_fbw=${wd}/results/tobias/ICM_ATAC_H3K27ac_filter3/TE_footprints.bw
# BINDetect
mkdir -p ${wd}/results/tobias/ICM_ATAC_H3K27ac_filter3/BIND_output_ICM_vs_TE
TOBIAS BINDetect --motifs ${jaspar} --signals ${icm_fbw} ${te_fbw} --genome ${genome_fa} --peaks ${binding_pks} \
                 --outdir ${wd}/results/tobias/ICM_ATAC_H3K27ac_filter3/BIND_output_ICM_vs_TE --cond_names ICM TE --cores 20
# PlotAggregate
mkdir -p ${wd}/results/tobias/ICM_ATAC_H3K27ac_filter3/PlotAggregate_output_ICM_vs_TE
ls ${wd}/results/tobias/ICM_ATAC_H3K27ac_filter3/BIND_output_ICM_vs_TE/ | while read dir
do
    tfbs_all=${wd}/results/tobias/ICM_ATAC_H3K27ac_filter3/BIND_output_ICM_vs_TE/${dir}/beds/${dir}_all.bed
    tfbs_icm_bound=${wd}/results/tobias/ICM_ATAC_H3K27ac_filter3/BIND_output_ICM_vs_TE/${dir}/beds/${dir}_ICM_bound.bed
    tfbs_icm_unbound=${wd}/results/tobias/ICM_ATAC_H3K27ac_filter3/BIND_output_ICM_vs_TE/${dir}/beds/${dir}_ICM_unbound.bed
    tfbs_te_bound=${wd}/results/tobias/ICM_ATAC_H3K27ac_filter3/BIND_output_ICM_vs_TE/${dir}/beds/${dir}_TE_bound.bed
    tfbs_te_unbound=${wd}/results/tobias/ICM_ATAC_H3K27ac_filter3/BIND_output_ICM_vs_TE/${dir}/beds/${dir}_TE_unbound.bed
    TOBIAS PlotAggregate --TFBS ${tfbs_all} ${tfbs_icm_bound} ${tfbs_icm_unbound} ${tfbs_te_bound} ${tfbs_te_unbound} \
                         --signals ${icm_cbw} ${naive_cbw1} ${naive_cbw2} ${primed_cbw1} ${primed_cbw2} ${te_cbw} \
                         --output ${wd}/results/tobias/ICM_ATAC_H3K27ac_filter3/PlotAggregate_output_ICM_vs_TE/${dir}_footprint_comparison_all.pdf \
                         --share_y both --plot_boundaries --signal-on-x --smooth 8 --flank 50
done



### >>> Find subfamily of repeat: open in ICM; close in TE; include motif of SOX2 and OCT4 (YHW)
motif=${wd}/metadata/motif.jaspar
repeat_consensus=/home/yhw/document/dfam/v3.3/hs_consensus/homo_sapiens_dfam_consensus.fa
npm2=./metadata/Npm2.fa
mkdir -p ${wd}/results/tobias/tfbscan_output
TOBIAS TFBScan --motifs ${motif} --fasta ${repeat_consensus} --cores 16 --outdir ${wd}/results/tobias/tfbscan_output
mkdir -p ${wd}/results/tobias/tfbscan_output_npm2
TOBIAS TFBScan --motifs ${motif} --fasta ./Npm2.fa --cores 16 --outdir ${wd}/results/tobias/tfbscan_output_npm2

