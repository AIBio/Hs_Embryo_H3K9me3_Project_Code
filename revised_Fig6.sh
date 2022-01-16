### >>> 1. Global setting
# - working directory
wd=/home/cmq/bioinfo/project-cmq/embryo_93
# - load fa file
grna_fa=${wd}/results/bedtools/gRNA_validate/sva_hervh.fa
genome_fa=/home/cmq/genome/ucsc/human/hg38/hg38.fa
# - load repeat annotation file
repeat_anno=/home/cmq/genome/ucsc/human/hg38/repeat/GRCh38_RepeatMasker_repeat_anno_quan.bed
# - load peak files
icm_43_pks=/home/cmq/bioinfo/project-cmq/gse124718_human_43/results/idr/broad0.01/SRR8402635_trimmed_psorted_dedup_peaks_SRR8402636_trimmed_psorted_dedup_peaks_ucsc_bl.bed
icm_93_pks=/home/cmq/bioinfo/project-cmq/embryo_93/results/idr/p_0.05/filtered/hICM93-IP-1_combined_R1_val_1_hICM93-IP-1_combined_R2_val_2_psorted_dedup_peaks_hICM93-IP-2_combined_R1_val_1_hICM93-IP-2_combined_R2_val_2_psorted_dedup_peaks_ucsc_bl.bed
primed_93_pks=/home/data/publicdata/GSE84382/analysis/results/idr/broad_p0.01/SRR3904792_trimmed_psort_dedup_peaks_SRR3904794_trimmed_psort_dedup_peaks.bed
naive_93_pks=/home/data/publicdata/GSE84382/analysis/results/idr/broad_p0.01/SRR3904795_trimmed_psort_dedup_peaks_SRR3904797_trimmed_psort_dedup_peaks.bed
icm_93_pks=/home/yhw/bioinfo/project-mine/Embryo.93/R/RawData/peaks/hs_k9/hICM93-IP-1_combined_R1_val_1_hICM93-IP-1_combined_R2_val_2_psorted_dedup_peaks_hICM93-IP-2_combined_R1_val_1_hICM93-IP-2_combined_R2_val_2_psorted_dedup_peaks_ucsc_bl.bed
naive_43_pks=/home/data/publicdata/GSE59434/analysis/results/macs2/narrow_p0.01/WIBR2_6iLA_H3K4me3_trimmed_psort_dedup_npks.bed
naive_273_pks=/home/data/publicdata/GSE59434/analysis/results/macs2/broad_p0.01/WIBR2_6iLA_H3K27me3_trimmed_psort_dedup_bpks.bed
primed_43_pks=/home/data/publicdata/GSE59434/analysis/results/macs2/narrow_p0.01/WIBR2_hESM_H3K4me3_trimmed_psort_dedup_npks.bed
primed_273_pks=/home/data/publicdata/GSE59434/analysis/results/macs2/broad_p0.01/WIBR2_hESM_H3K27me3_trimmed_psort_dedup_bpks.bed
icm_273_pks=/home/cmq/bioinfo/project-cmq/gse124718/results/idr/multi_merged_BAMPE_p_0.01/SRR8402638_trimmed_psorted_dedup_peaks_SRR8402639_trimmed_psorted_dedup_peaks_ucsc_bl.bed
icm_43_pks=/home/cmq/bioinfo/project-cmq/gse124718_human_43/results/idr/broad0.01/SRR8402635_trimmed_psorted_dedup_peaks_SRR8402636_trimmed_psorted_dedup_peaks_ucsc_bl.bed
pro1k=/home/yhw/genome/ensembl/release97/homo_sapiens/dna_anno/regions/Homo_sapiens.GRCh38.97_pro1k.bed
hs_repeat=./metadata/GRCh38_RepeatMasker_repeat_region.bed
hs_repeat_plus_pro1k=./metadata/GRCh38_RepeatMasker_repeat_and_gene_pro1k.bed
cat ${icm_93_pks} ${naive_93_pks} ${primed_93_pks} | sort -k1,1 -k2,2n | cut -f 1-3 | bedtools merge -i stdin > ./results/bedtools/naive/ICM_Naive_Primed_merge_peaks.bed
both_merged_93pks=./results/bedtools/naive/ICM_Naive_Primed_merge_peaks.bed
cat ${icm_93_pks} ${naive_93_pks} ${primed_93_pks} | sort -k1,1 -k2,2n | cut -f 1-3 | bedtools merge -i stdin > ./results/bedtools/naive/ICM_Naive_Primed_merge_peaks.bed
both_merged_93pks=./results/bedtools/naive/ICM_Naive_Primed_merge_peaks.bed
hs_sva=/home/data/wanglab/ICM_like/analysis/chipseq/metadata/GRCh38_RepeatMasker_all_sva.bed

# - load bw files
primed_93_input=/home/data/publicdata/GSE84382/analysis/results/deeptools/coverage/SRR3904789_trimmed_psort_dedup.bw
naive_93_input=/home/data/publicdata/GSE84382/analysis/results/deeptools/coverage/SRR3904791_trimmed_psort_dedup.bw
primed_93_1=/home/data/publicdata/GSE84382/analysis/results/deeptools/coverage/SRR3904792_trimmed_psort_dedup.bw
primed_93_2=/home/data/publicdata/GSE84382/analysis/results/deeptools/coverage/SRR3904794_trimmed_psort_dedup.bw
naive_93_1=/home/data/publicdata/GSE84382/analysis/results/deeptools/coverage/SRR3904795_trimmed_psort_dedup.bw
naive_93_2=/home/data/publicdata/GSE84382/analysis/results/deeptools/coverage/SRR3904797_trimmed_psort_dedup.bw

naive_273=/home/data/publicdata/GSE59434/analysis/results/deeptools/coverage/dedup_bam/WIBR2_6iLA_H3K27me3_trimmed_psort_dedup.bw
naive_43=/home/data/publicdata/GSE59434/analysis/results/deeptools/coverage/dedup_bam/WIBR2_6iLA_H3K4me3_trimmed_psort_dedup.bw
naive_43=/home/data/publicdata/GSE59434/analysis/results/deeptools/coverage/dedup_bam_scale2/WIBR2_6iLA_H3K4me3_trimmed_psort_dedup.bw
naive_43_273_input=/home/data/publicdata/GSE59434/analysis/results/deeptools/coverage/dedup_bam/WIBR2_6iLA_Input_trimmed_psort_dedup.bw
naive_273_fc=/home/data/publicdata/GSE59434/analysis/results/deeptools/bamcompare/dedup_bam/WIBR2_6iLA_H3K27me3_trimmed_psort_dedup_vs_WIBR2_6iLA_Input_trimmed_psort_dedup.bw
naive_43_fc=/home/data/publicdata/GSE59434/analysis/results/deeptools/bamcompare/dedup_bam/WIBR2_6iLA_H3K4me3_trimmed_psort_dedup_vs_WIBR2_6iLA_Input_trimmed_psort_dedup.bw

primed_273=/home/data/publicdata/GSE59434/analysis/results/deeptools/coverage/dedup_bam/WIBR2_hESM_H3K27me3_trimmed_psort_dedup.bw
primed_43=/home/data/publicdata/GSE59434/analysis/results/deeptools/coverage/dedup_bam/WIBR2_hESM_H3K4me3_trimmed_psort_dedup.bw
primed_43=/home/data/publicdata/GSE59434/analysis/results/deeptools/coverage/dedup_bam_scale2/WIBR2_hESM_H3K4me3_trimmed_psort_dedup.bw
primed_43_273_input=/home/data/publicdata/GSE59434/analysis/results/deeptools/coverage/dedup_bam/WIBR2_hESM_Input_trimmed_psort_dedup.bw
primed_273_fc=/home/data/publicdata/GSE59434/analysis/results/deeptools/bamcompare/dedup_bam/WIBR2_hESM_H3K27me3_trimmed_psort_dedup_vs_WIBR2_hESM_Input_trimmed_psort_dedup.bw
primed_43_fc=/home/data/publicdata/GSE59434/analysis/results/deeptools/bamcompare/dedup_bam/WIBR2_hESM_H3K4me3_trimmed_psort_dedup_vs_WIBR2_hESM_Input_trimmed_psort_dedup.bw

primed_93_fc1=/home/data/publicdata/GSE84382/analysis/results/deeptools/coverage/foldchange/SRR3904792_trimmed_psort_dedup_vs_SRR3904789_trimmed_psort_dedup.bw
primed_93_fc2=/home/data/publicdata/GSE84382/analysis/results/deeptools/coverage/foldchange/SRR3904794_trimmed_psort_dedup_vs_SRR3904789_trimmed_psort_dedup.bw
naive_93_fc1=/home/data/publicdata/GSE84382/analysis/results/deeptools/coverage/foldchange/SRR3904795_trimmed_psort_dedup_vs_SRR3904791_trimmed_psort_dedup.bw
naive_93_fc2=/home/data/publicdata/GSE84382/analysis/results/deeptools/coverage/foldchange/SRR3904797_trimmed_psort_dedup_vs_SRR3904791_trimmed_psort_dedup.bw

icm_input=/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/hs_k9me3/results/deeptools/coverage/merge/hICM93-Input.bw
icm_ip=/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/hs_k9me3/results/deeptools/coverage/merge/hICM93-IP-1.bw
icm_fc=/home/cmq/bioinfo/project-cmq/embryo_93/results/deeptools/bigwig/merged_bam_ratio/hICM93-IP-1_vs_hICM93-Input.bw
icm_k273=/home/cmq/bioinfo/project-cmq/gse124718/results/deeptools/bigwig/merge_replicate_multi/icm_multi_replicate_merged.bw
icm_k43=/home/cmq/bioinfo/project-cmq/gse124718_human_43/results/deeptools/bigwig/merged_bam/SRR8402635_trimmed_psorted_dedup_SRR8402636_trimmed_psorted_dedup_merged.bw

ln -s ${icm_ip} ./results/deeptools/coverage/naive/ICM_93.bw
ln -s ${icm_fc} ./results/deeptools/coverage/naive/ICM_93_lfc.bw
ln -s ${naive_93_fc1} ./results/deeptools/coverage/naive/Naive_93_lfc.bw
ln -s ${naive_93_1} ./results/deeptools/coverage/naive/Naive_93.bw
ln -s ${primed_93_fc1} ./results/deeptools/coverage/naive/Primed_93_lfc.bw
ln -s ${primed_93_1} ./results/deeptools/coverage/naive/Primed_93.bw
ln -s ${icm_k273} ./results/deeptools/coverage/naive/ICM_273.bw
ln -s ${icm_k43} ./results/deeptools/coverage/naive/ICM_43.bw
ln -s ${naive_273} ./results/deeptools/coverage/naive/Naive_273.bw
ln -s ${naive_43} ./results/deeptools/coverage/naive/Naive_43.bw
ln -s ${primed_273} ./results/deeptools/coverage/naive/Primed_273.bw
ln -s ${primed_43} ./results/deeptools/coverage/naive/Primed_43.bw

c4_input=results/deeptools/coverage/merge/h4C93-Input-1.bw
c4_ip=results/deeptools/coverage/merge/h4C93-Ip-1.bw
c4_fc=/home/cmq/bioinfo/project-cmq/embryo_93/results/deeptools/bigwig/merged_bam_ratio/h4C93-Ip-1_vs_h4C93-Input-1.bw

c8_input=results/deeptools/coverage/merge/h8C93-Input-2.bw
c8_ip=results/deeptools/coverage/merge/h8C93-Ip-3.bw
c8_fc=/home/cmq/bioinfo/project-cmq/embryo_93/results/deeptools/bigwig/merged_bam_ratio/h8C93-Ip-3_vs_h8C93-Input-2.bw

morula_input=results/deeptools/coverage/merge/hM93-Input.bw
morula_ip=results/deeptools/coverage/merge/hM93-Ip-1.bw
morula_fc=/home/cmq/bioinfo/project-cmq/embryo_93/results/deeptools/bigwig/merged_bam_ratio/hM93-Ip-1_vs_hM93-Input.bw

icm_input=results/deeptools/coverage/merge/hICM93-Input.bw
icm_ip=results/deeptools/coverage/merge/hICM93-IP-1.bw
icm_fc=/home/cmq/bioinfo/project-cmq/embryo_93/results/deeptools/bigwig/merged_bam_ratio/hICM93-IP-1_vs_hICM93-Input.bw

c4_43_1=/home/cmq/bioinfo/project-cmq/gse124718_human_43/results/deeptools/bigwig/merged_bam/SRR8402611_1_val_1_SRR8402611_2_val_2_psorted_dedup_SRR8402612_trimmed_psorted_dedup_merged.bw
c4_43_2=/home/cmq/bioinfo/project-cmq/gse124718_human_43/results/deeptools/bigwig/merged_bam/SRR8402613_trimmed_psorted_dedup_SRR8402614_1_val_1_SRR8402614_2_val_2_psorted_dedup_merged.bw
c8_43=/home/cmq/bioinfo/project-cmq/gse124718_human_43/results/deeptools/bigwig/merged_bam/SRR8402624_1_val_SRR8402625_1_val_SRR8402627_1_val_SRR8402628_1_valmerged.bw
icm_43=/home/cmq/bioinfo/project-cmq/gse124718_human_43/results/deeptools/bigwig/merged_bam/SRR8402635_trimmed_psorted_dedup_SRR8402636_trimmed_psorted_dedup_merged.bw

c2_273=/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/hs_k27me3/results/deeptools/coverage/merge/2cell_273.bw
c4_273=/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/hs_k27me3/results/deeptools/coverage/merge/4cell_273.bw
c8_273=/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/hs_k27me3/results/deeptools/coverage/merge/8cell_273.bw
morula_273=/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/hs_k27me3/results/deeptools/coverage/merge/morula_273.bw
icm_273=/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/hs_k27me3/results/deeptools/coverage/merge/icm_273.bw
te_273=/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/hs_k27me3/results/deeptools/coverage/merge/te_273.bw

c8_27ac=/home/yhw/bioinfo/project-temp/wjc/gse124718/results/deeptools/coverage/SRR9131735_1_val_1_SRR9131735_2_val_2_psort_dedup.bw
icm_27ac=/home/yhw/bioinfo/project-temp/wjc/gse124718/results/deeptools/coverage/SRR9131737_1_val_1_SRR9131737_2_val_2_psort_dedup.bw
naive_27ac=/home/data/publicdata/GSE69646/analysis/results/deeptools/coverage/dedup_bam/Naive_H3K27ac_rep1_trimmed_psort_dedup.bw
primed_27ac=/home/data/publicdata/GSE69646/analysis/results/deeptools/coverage/dedup_bam/Primed_H3K27ac_rep1_trimmed_psort_dedup.bw



### >>> 2. Marked gene and repeat: ICM, Naive and Primed H3K9me3 H3K27me3 H3K4me3 (CMQ)
mkdir -p results/bedtools/naive/marked_repeat_gene
for pks in icm_93_pks icm_273_pks icm_43_pks naive_93_pks naive_273_pks naive_43_pks primed_93_pks primed_273_pks primed_43_pks
do
    wc -l ${!pks}
    bedtools intersect -a ${hs_repeat_plus_pro1k} -b ${!pks} -wa -f 0.5 > ./results/bedtools/naive/marked_repeat_gene/${pks}_marked_repeat_and_pro1k.bed
done



### >>> 3. Validation of SVA gRNA (CMQ)
# - mapping gRNA to reference genome
od=${wd}/results/bedtools/gRNA_validate
pblat -t=dna -q=dna -minScore=0 -stepSize=2 -oneOff=4 -threads=12 ${genome_fa} ${grna_fa} ${od}/gRNA_validate.psl
# - extract the potentially targeted regions
cat ${od}/gRNA_validate.psl | awk '{if($2<=2 && $7==0 && ($13-$12)>=15) print $14"\t"$16"\t"$17"\t"$10}' > ${od}/gRNA_target_region.bed
# - extract the targeted repeats
grep "SVA" ${repeat_anno} > ${od}/all_SVA.bed
bedtools intersect -a ${od}/gRNA_target_region.bed -b ${repeat_anno} -f 1 -wa -wb | grep ":SVA:" > ${od}/gRNA_target_SVA.bed
bedtools intersect -a ${od}/gRNA_target_region.bed -b ${repeat_anno} -f 1 -wa -wb > ${od}/gRNA_target_repeat.bed



### >>> 4. Bivalent SVAs (K9&K4) (YHW)
# - define bivalent peaks
icm_43_pks=/home/cmq/bioinfo/project-cmq/gse124718_human_43/results/idr/broad0.05/SRR8402635_trimmed_psorted_dedup_peaks_SRR8402636_trimmed_psorted_dedup_peaks_ucsc_bl.bed
bedtools intersect -a ${icm_93_pks} -b ${icm_43_pks} \
   | awk '{if(($3-$2)>100) print $0}' \
   | bedtools intersect -a stdin -b ${icm_vs_naive_down_repeat} -wa > ./results/bedtools/naive/icm_vs_naive_down_repeat_marked_by_bivalent_K93_K43.bed
bedtools intersect -a ${icm_93_pks} -b ${icm_43_pks} \
   | awk '{if(($3-$2)>100) print $0}' \
   | bedtools intersect -a stdin -b ${icm_vs_naive_up_repeat} -wa > ./results/bedtools/naive/icm_vs_naive_up_repeat_marked_by_bivalent_K93_K43.bed
# - define H3K9me3 and H3K4me3 overlapped regions
od=${wd}/results/bedtools/revised_k9k4_domains
bedtools intersect -a ${icm_93_pks} -b ${icm_43_pks} | awk '{if(($3-$2)>200) print $0}' > ${od}/hICM_H3K9me3_and_H3K4me3_bivalent_peaks.bed
# - extract H3K9me3 and H3K4me3 overlapped SVA
bedtools intersect -a ${repeat_anno} -b ${od}/hICM_H3K9me3_and_H3K4me3_bivalent_peaks.bed -f 0.5 -wa \
   | grep ":SVA:" | grep -v "fix\|alt" | sort -k1,1 -k2,2n > ${od}/hICM_H3K9me3_and_H3K4me3_bivalent_peaks_marked_SVA_0.5pect.bed
# - plotting heatmap
for pks in icm_bi_SVA icm_bi_SVA_LTR icm_bi_SVA_LTR_L1
do
    mkdir -p ./results/deeptools/matrix/scale5000/naive ./results/deeptools/heatmap/scale5000/naive
    computeMatrix scale-regions -S ${c4_fc} ${c4_43_1} ${c4_43_2} ${c4_273} ${c8_fc} ${c8_43} ${c8_273} ${c8_27ac} ${icm_fc} ${icm_43} ${icm_273} ${icm_27ac} ${naive_93_fc1} ${naive_43} ${naive_273} ${naive_27ac} ${primed_93_fc1} ${primed_43} ${primed_273} ${primed_27ac} \
                                -R ${!pks} \
                                --binSize 20 --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
                                --numberOfProcessors 40 --missingDataAsZero -o ./results/deeptools/matrix/scale5000/naive/${pks}_RPKM_and_LFC.matrix.gz
    plotHeatmap -m ./results/deeptools/matrix/scale5000/naive/${pks}_RPKM_and_LFC.matrix.gz \
                -o ./results/deeptools/heatmap/scale5000/naive/${pks}_K9_LFC.matrix.pdf \
                --averageTypeSummaryPlot mean --colorList "cyan, white, red" --heatmapHeight 6 --heatmapWidth 2 --whatToShow "plot, heatmap and colorbar" \
                --plotType lines --zMin -1 --zMax 5 --yMin -1 --yMax 3 --legendLocation upper-left --interpolationMethod auto \
                --samplesLabel "4C 93" "4C 43-1" "4C 43-2" "4C 273" "8C 93" "8C 43" "8C 273" "8C 27ac" "ICM 93" "ICM 43" "ICM 273" "ICM 27ac" "Naive 93" "Naive 43" "Naive 273" "Naive 27ac" "Primed 93" "Primed 43" "Primed 273" "Primed 27ac"
    plotHeatmap -m ./results/deeptools/matrix/scale5000/naive/${pks}_RPKM_and_LFC.matrix.gz \
                -o ./results/deeptools/heatmap/scale5000/naive/${pks}_K4_embryo_RPKM.matrix.pdf \
                --averageTypeSummaryPlot mean --colorMap Blues --heatmapHeight 6 --heatmapWidth 2 --whatToShow "plot, heatmap and colorbar" \
                --plotType lines --zMin 0 --zMax 100 --legendLocation upper-left --interpolationMethod auto \
                --samplesLabel "4C 93" "4C 43-1" "4C 43-2" "4C 273" "8C 93" "8C 43" "8C 273" "8C 27ac" "ICM 93" "ICM 43" "ICM 273" "ICM 27ac" "Naive 93" "Naive 43" "Naive 273" "Naive 27ac" "Primed 93" "Primed 43" "Primed 273" "Primed 27ac"
    plotHeatmap -m ./results/deeptools/matrix/scale5000/naive/${pks}_RPKM_and_LFC.matrix.gz \
                -o ./results/deeptools/heatmap/scale5000/naive/${pks}_K4_naive_primed_RPKM.matrix.pdf \
                --averageTypeSummaryPlot mean --colorMap Blues --heatmapHeight 6 --heatmapWidth 2 --whatToShow "plot, heatmap and colorbar" \
                --plotType lines --zMin 0 --zMax 50 --yMin 0 --yMax 50 --legendLocation upper-left --interpolationMethod auto \
                --samplesLabel "4C 93" "4C 43-1" "4C 43-2" "4C 273" "8C 93" "8C 43" "8C 273" "8C 27ac" "ICM 93" "ICM 43" "ICM 273" "ICM 27ac" "Naive 93" "Naive 43" "Naive 273" "Naive 27ac" "Primed 93" "Primed 43" "Primed 273" "Primed 27ac"
    plotHeatmap -m ./results/deeptools/matrix/scale5000/naive/${pks}_RPKM_and_LFC.matrix.gz \
                -o ./results/deeptools/heatmap/scale5000/naive/${pks}_K27_embryo_RPKM.matrix.pdf \
                --averageTypeSummaryPlot mean --colorMap Blues --heatmapHeight 6 --heatmapWidth 2 --whatToShow "plot, heatmap and colorbar" \
                --plotType lines --zMin 0 --zMax 30 --yMin 0 --yMax 30 --legendLocation upper-left --interpolationMethod auto \
                --samplesLabel "4C 93" "4C 43-1" "4C 43-2" "4C 273" "8C 93" "8C 43" "8C 273" "8C 27ac" "ICM 93" "ICM 43" "ICM 273" "ICM 27ac" "Naive 93" "Naive 43" "Naive 273" "Naive 27ac" "Primed 93" "Primed 43" "Primed 273" "Primed 27ac"
    plotHeatmap -m ./results/deeptools/matrix/scale5000/naive/${pks}_RPKM_and_LFC.matrix.gz \
                -o ./results/deeptools/heatmap/scale5000/naive/${pks}_K27_naive_primed_RPKM.matrix.pdf \
                --averageTypeSummaryPlot mean --colorMap Blues --heatmapHeight 6 --heatmapWidth 2 --whatToShow "plot, heatmap and colorbar" \
                --plotType lines --zMin 0 --zMax 5 --yMin 0 --yMax 5 --legendLocation upper-left --interpolationMethod auto \
                --samplesLabel "4C 93" "4C 43-1" "4C 43-2" "4C 273" "8C 93" "8C 43" "8C 273" "8C 27ac" "ICM 93" "ICM 43" "ICM 273" "ICM 27ac" "Naive 93" "Naive 43" "Naive 273" "Naive 27ac" "Primed 93" "Primed 43" "Primed 273" "Primed 27ac"
    plotHeatmap -m ./results/deeptools/matrix/scale5000/naive/${pks}_RPKM_and_LFC.matrix.gz \
                -o ./results/deeptools/heatmap/scale5000/naive/${pks}_K27ac_embryo_RPKM.matrix.pdf \
                --averageTypeSummaryPlot mean --colorMap Blues --heatmapHeight 6 --heatmapWidth 2 --whatToShow "plot, heatmap and colorbar" \
                --plotType lines --zMin 0 --zMax 40 --yMin 0 --yMax 30 --legendLocation upper-left --interpolationMethod auto \
                --samplesLabel "4C 93" "4C 43-1" "4C 43-2" "4C 273" "8C 93" "8C 43" "8C 273" "8C 27ac" "ICM 93" "ICM 43" "ICM 273" "ICM 27ac" "Naive 93" "Naive 43" "Naive 273" "Naive 27ac" "Primed 93" "Primed 43" "Primed 273" "Primed 27ac"
done



### >>> 5. cell line specific H3K9me3 peaks (YHW)
HomerMotif(){        
   indir=$1; file=$2; outdir=$3; logdir=$4; gv=$5
   cat ${file} | while read region
   do
       prefix="${region%.*}"
       findMotifsGenome.pl ${indir}/${region} ${gv} ${outdir} > ${logdir}/${prefix}.log 2>&1
   done   
}
# - ICM vs Naive
icm_vs_naive_93_pks=./results/bedtools/naive/ICM_vs_Naive_specific_pks_p0.05.bed
bedtools intersect -a ${icm_93_pks} -b ${naive_93_pks} -v -wa | sort | uniq | sort -k1,1 -k2,2n > ./results/bedtools/naive/ICM_vs_Naive_specific_pks_p0.05.bed
bedtools intersect -a metadata/GRCh38_RepeatMasker_repeat_and_gene_region.bed -b ${icm_vs_naive_93_pks} -wa -f 0.5 > ./results/bedtools/naive/ICM_vs_Naive_specific_pks_p0.05_marked_gene_repeat.bed
# - ICM vs Primed
icm_vs_primed_93_pks=./results/bedtools/naive/ICM_vs_Primed_specific_pks_p0.05.bed
bedtools intersect -a ${icm_93_pks} -b ${primed_93_pks} -v -wa | sort | uniq | sort -k1,1 -k2,2n > ./results/bedtools/naive/ICM_vs_Primed_specific_pks_p0.05.bed
# - Naive vs ICM
naive_vs_icm_93_pks=./results/bedtools/naive/Naive_vs_ICM_specific_pks_p0.05.bed
bedtools intersect -a ${naive_93_pks} -b ${icm_93_pks} -v -wa | sort | uniq | sort -k1,1 -k2,2n > ./results/bedtools/naive/Naive_vs_ICM_specific_pks_p0.05.bed
bedtools intersect -a metadata/GRCh38_RepeatMasker_repeat_and_gene_region.bed -b ${naive_vs_icm_93_pks} -wa -f 0.5 > ./results/bedtools/naive/Naive_vs_ICM_specific_pks_p0.05_marked_gene_repeat.bed
# - Naive vs Primed
naive_vs_primed_93_pks=./results/bedtools/naive/Naive_vs_Primed_specific_pks_p0.05.bed
bedtools intersect -a ${naive_93_pks} -b ${primed_93_pks} -v -wa | sort | uniq | sort -k1,1 -k2,2n > ./results/bedtools/naive/Naive_vs_Primed_specific_pks_p0.05.bed
# - Primed vs ICM
primed_vs_icm_93_pks=./results/bedtools/naive/Primed_vs_ICM_specific_pks_p0.05.bed
bedtools intersect -a ${primed_93_pks} -b ${icm_93_pks} -v -wa | sort | uniq | sort -k1,1 -k2,2n > ./results/bedtools/naive/Primed_vs_ICM_specific_pks_p0.05.bed
# - Primed vs Naive
primed_vs_naive_93_pks=./results/bedtools/naive/Primed_vs_Naive_specific_pks_p0.05.bed
bedtools intersect -a ${primed_93_pks} -b ${naive_93_pks} -v -wa | sort | uniq | sort -k1,1 -k2,2n > ./results/bedtools/naive/Primed_vs_Naive_specific_pks_p0.05.bed
# - ICM, Naive and Primed shared
both_shared_93_shared=./results/bedtools/naive/ICM_Naive_Primed_shared_pks_p0.05.bed
bedtools intersect -a ${icm_93_pks} -b ${naive_93_pks} \
   | bedtools intersect -a stdin -b ${primed_93_pks} \
   | awk '{if(($3-$2)>300) print $0}' | sort -k1,1 -k2,2n > ./results/bedtools/naive/ICM_Naive_Primed_shared_pks_p0.05.bed
# - Plot heatmap
mkdir -p ./results/deeptools/matrix/reference/naive ./results/deeptools/heatmap/reference/naive
computeMatrix reference-point -S ${icm_fc} ${naive_93_fc1} ${primed_93_fc1} \
                              -R ${icm_vs_naive_93_pks} ${icm_vs_primed_93_pks} ${naive_vs_icm_93_pks} ${naive_vs_primed_93_pks} ${primed_vs_icm_93_pks} ${primed_vs_naive_93_pks} \
                              --averageTypeBins mean --referencePoint center \
                              --binSize 100 --beforeRegionStartLength 3000 --afterRegionStartLength 3000 \
                              --numberOfProcessors 30 --missingDataAsZero -o ./results/deeptools/matrix/reference/naive/icm_naive_primed_spe_93_pks.matrix.gz
plotHeatmap -m ./results/deeptools/matrix/reference/naive/icm_naive_primed_spe_93_pks.matrix.gz \
            -o ./results/deeptools/heatmap/reference/naive/icm_naive_primed_spe_93_pks.matrix.pdf \
            --averageTypeSummaryPlot mean --colorList "cyan, white, red" --heatmapHeight 8 --heatmapWidth 2 --whatToShow "plot, heatmap and colorbar" \
            --plotType lines --zMin -1 --zMax 3 --legendLocation upper-left --interpolationMethod auto \
            --samplesLabel "ICM" "Naive" "Primed"
for pks in naive_93_pks primed_93_pks icm_93_pks
do
    computeMatrix reference-point -S ${icm_fc} ${naive_93_fc1} ${primed_93_fc1} \
                                  -R ${!pks} \
                                  --averageTypeBins mean --referencePoint center \
                                  --binSize 100 --beforeRegionStartLength 3000 --afterRegionStartLength 3000 \
                                  --numberOfProcessors 30 --missingDataAsZero -o ./results/deeptools/matrix/reference/naive/${pks}.matrix.gz
    plotHeatmap -m ./results/deeptools/matrix/reference/naive/${pks}.matrix.gz \
                -o ./results/deeptools/heatmap/reference/naive/${pks}.matrix.pdf \
                --averageTypeSummaryPlot mean --colorList "cyan, white, red" --heatmapHeight 6 --heatmapWidth 2 --whatToShow "plot, heatmap and colorbar" \
                --plotType lines --zMin -1 --zMax 3 --legendLocation upper-left --interpolationMethod auto \
                --samplesLabel "ICM" "Naive" "Primed"
done



### >>> 6. Plotting H3K4me3 peaks (YHW)
for pks in naive_273_pks primed_273_pks icm_43_pks icm_273_pks naive_43_pks primed_43_pks 
do
    # ICM: H3K27me3 & H3K4me3 RPKM
    # ICM: H3K9me3 lfc; Naive & Primed H3K9me3 & H3K27me3 & H3K4me3 lfc
    computeMatrix reference-point -S ${icm_fc} ${naive_93_fc1} ${primed_93_fc1} \
                                     ${icm_k273} ${naive_273} ${naive_273_fc} ${primed_273} ${primed_273_fc} \
                                     ${icm_k43} ${naive_43} ${naive_43_fc} ${primed_43} ${primed_43_fc} \
                                  -R ${!pks} \
                                  --averageTypeBins mean --referencePoint center \
                                  --binSize 100 --beforeRegionStartLength 3000 --afterRegionStartLength 3000 \
                                  --numberOfProcessors 40 --missingDataAsZero -o ./results/deeptools/matrix/reference/naive/${pks}_RPKM_and_LFC.matrix.gz
    plotHeatmap -m ./results/deeptools/matrix/reference/naive/${pks}_RPKM_and_LFC.matrix.gz \
                -o ./results/deeptools/heatmap/reference/naive/${pks}_RPKM.matrix.pdf \
                --averageTypeSummaryPlot mean --colorMap Blues --heatmapHeight 6 --heatmapWidth 2 --whatToShow "plot, heatmap and colorbar" \
                --plotType lines --zMin 0 --zMax 20 --legendLocation upper-left --interpolationMethod auto \
                --samplesLabel "ICM 93 FC" "Naive 93 FC" "Primed 93 FC" \
                               "ICM 273 R" "Naive 273 R" "Naive 273 FC" "Primed 273 R" "Primed 273 FC" \
                               "ICM 43" "Naive 43 R" "Naive 43 FC" "Primed 43 R" "Primed 43 FC"
    plotHeatmap -m ./results/deeptools/matrix/reference/naive/${pks}_RPKM_and_LFC.matrix.gz \
                -o ./results/deeptools/heatmap/reference/naive/${pks}_LFC.matrix.pdf \
                --averageTypeSummaryPlot mean --colorList "cyan, white, red" --heatmapHeight 6 --heatmapWidth 2 --whatToShow "plot, heatmap and colorbar" \
                --plotType lines --zMin -1 --zMax 2 --yMin 0 --yMax 3 --legendLocation upper-left --interpolationMethod auto \
                --samplesLabel "ICM 93 FC" "Naive 93 FC" "Primed 93 FC" \
                               "ICM 273 R" "Naive 273 R" "Naive 273 FC" "Primed 273 R" "Primed 273 FC" \
                               "ICM 43" "Naive 43 R" "Naive 43 FC" "Primed 43 R" "Primed 43 FC"
done



### >>> 7. Calculate values in specific regions (YHW)
RegionCoverageNor(){
   indir=$1; outdir=$2; logdir=$3; key=$4; region=$5
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   bwlist=$(ls ${indir}/*${key}*.bw | tr '\n' ' ')
   fileName="${region##*/}"; prefix="${fileName%.*}"
   multiBigwigSummary BED-file --BED ${region} -p 16 -b ${bwlist} -o ${outdir}/${prefix}_readNor.npz \
                               --outRawCounts ${outdir}/${prefix}_readNor.tab > ${logdir}/${prefix}_nor.log 2>&1
}
RegionCoverageNor ./results/deeptools/coverage/naive ./results/deeptools/multibw/naive_icm_93pks ./results/deeptools/multibw/naive_icm_93pks "*" ${icm_93_pks}
RegionCoverageNor ./results/deeptools/coverage/naive ./results/deeptools/multibw/naive_icm_pro1k ./results/deeptools/multibw/naive_icm_pro1k "*" ${pro1k}
RegionCoverageNor ./results/deeptools/coverage/naive ./results/deeptools/multibw/naive_icm_93pks ./results/deeptools/multibw/naive_icm_93pks "*" ${both_merged_93pks}
RegionCoverageNor ./results/deeptools/coverage/naive ./results/deeptools/multibw/naive_icm_repeat ./results/deeptools/multibw/naive_icm_repeat "*" ${hs_repeat}
RegionCoverageNor ./results/deeptools/coverage/naive ./results/deeptools/multibw/naive_icm_hs_repeat_plus_pro1k ./results/deeptools/multibw/naive_icm_hs_repeat_plus_pro1k "*" ${hs_repeat_plus_pro1k}
RegionCoverageNor ./results/deeptools/coverage/naive ./results/deeptools/multibw/naive_icm_sva ./results/deeptools/multibw/naive_icm_sva "*" ${hs_sva}
GlobalCoverageNor(){
   indir=$1; outdir=$2; logdir=$3; key=$4; prefix=$5
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   bwlist=$(ls ${indir}/*${key}*.bw | tr '\n' ' ')
   multiBigwigSummary bins -bs 10000 -n 0 -p 16 -b ${bwlist} -o ${outdir}/${prefix}_readNor.npz \
                           --outRawCounts ${outdir}/${prefix}_readNor.tab > ${logdir}/${prefix}_nor.log 2>&1
}
GlobalCoverageNor ./results/deeptools/coverage/naive ./results/deeptools/multibw/naive_icm_global ./logs/deeptools/multibw/naive_icm_global "*" "global_signal"

