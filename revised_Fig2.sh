### >>> Objects
# - Global setting
# - Define enhancers: RE.by.H3K9me3, RE.by.H3K27me3 and PE(primed) (YHW)
# - Plotting signal on PE and RE (YHW)
# - Footprint analysis (YHW)
# - Find the genes interacted by self-defined enhancers: RE.K9, RE.K273, PE (YHW)
# - Find the genes interacted by self-defined enhancers: RE.K9, RE.K273, PE (YHW)
# - Interaction aggregate plots (YHW)
# - Plotting signal on mouse ZGA genes (CMQ)
# - Plotting signal on human ZGA genes (CMQ)
# - Plot Hic aggregate (YHW)
# - Hic data contact frequency with ZGA genes (YHW)



### >>> Global setting
# - load peaks
pro3k=/home/cmq/genome/ensembl/release97/homo_sapiens/Homo_sapiens.GRCh38.97.ens2ucsc.all.gene.pro3k.bed
pro1k=/home/cmq/genome/ensembl/release97/homo_sapiens/Homo_sapiens.GRCh38.97.ens2ucsc.all.gene.pro1k.bed
c8_27ac_pks=/home/yhw/bioinfo/project-temp/wjc/gse124718/results/idr/narrow_p0.05/SRR9131735_1_val_1_SRR9131735_2_val_2_psort_dedup_peaks_SRR9131736_1_val_1_SRR9131736_2_val_2_psort_dedup_peaks.bed
c8_atac_pks=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/idr/dedup_bam/p_0.05/SRR5837329_SRR5837329_psorted_dedup_SRR5837330_SRR5837330_psorted_dedup_merged_peaks_SRR5837331_SRR5837331_psorted_dedup_SRR5837332_SRR5837332_psorted_dedup_merged_peaks_ucsc_bl.bed
c8_loops=/home/yhw/bioinfo/project-mine/Embryo.93/hic/hs_embryo/results/fithic/8cell_rep2_loops/FitHiC.spline_pass1.res40000.significances.txt.gz
c4_93_pks=/home/cmq/bioinfo/project-cmq/embryo_93/results/idr/p_0.05/filtered/h4C93-Ip-1_combined_R1_val_1_h4C93-Ip-1_combined_R2_val_2_psorted_dedup_peaks_h4C93-Ip-2_combined_R1_val_1_h4C93-Ip-2_combined_R2_val_2_psorted_dedup_peaks_ucsc_bl.bed
c4_atac_pks=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/idr/dedup_bam/p_0.05/SRR5837324_SRR5837324_psorted_dedup_SRR5837325_SRR5837325_psorted_dedup_merged_peaks_SRR5837326_SRR5837326_psorted_dedup_SRR5837327_SRR5837327_psorted_dedup_SRR5837328_SRR5837328_psorted_dedup_merged_peaks_ucsc_bl.bed
c2_273_pks=/home/cmq/bioinfo/project-cmq/gse124718/results/idr/multi_merged_BAMPE_p_0.05/filtered/SRR8402607_SRR8402607_psorted_dedup_SRR8402608_SRR8402608_psorted_dedup_merged_peaks_SRR8402609_SRR8402609_psorted_dedup_SRR8402610_SRR8402610_psorted_dedup_merged_peaks_ucsc_bl.bed
c4_273_pks=/home/cmq/bioinfo/project-cmq/gse124718/results/idr/multi_merged_BAMPE_p_0.05/filtered/SRR8402616_SRR8402616_psorted_dedup_peaks_SRR8402617_SRR8402617_psorted_dedup_peaks_ucsc_bl.bed
c8_273_pks=/home/cmq/bioinfo/project-cmq/gse124718/results/idr/multi_merged_BAMPE_p_0.05/filtered/SRR8402630_SRR8402630_psorted_dedup_SRR8402631_SRR8402631_psorted_dedup_merged_peaks_SRR8402633_SRR8402633_psorted_dedup_SRR8402634_SRR8402634_psorted_dedup_merged_peaks_ucsc_bl.bed



### >>> Define enhancers: RE.by.H3K9me3, RE.by.H3K27me3 and PE(primed) (YHW)
od=/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer
[ ! -e "${od}" ] && mkdir -p ${od}
# - ChIPseq&ATACseq: extract distal H3K27ac and ATAC peaks
bedtools intersect -a ${c8_27ac_pks} -b ${pro3k} -v | awk '{if(($3-$2)>=200) print $0}' | cut -f 1-3 > ${od}/8cell_distal_h3k27ac_pks.bed
bedtools intersect -a ${c8_atac_pks} -b ${pro3k} -v | awk '{if(($3-$2)>=200) print $0}' | cut -f 1-3 > ${od}/8cell_distal_atac_pks.bed

# - HiC: extract loops contained genes
zcat ${c8_loops} | awk '{OFS=FS="\t"}{if($6<=0.05&&$5>=20) print $1,$2-20000,$2+20000,$3,$4-20000,$4+20000,$5,$6}' | sort -k1,1 -k2,2n > ${od}/8cell_significant_loops.bed
loops_bin=${od}/8cell_significant_loops_bin.bed
awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3"\n"$4,$5,$6}' ${od}/8cell_significant_loops.bed | sort -u > ${loops_bin}

# - HiC: extract loops overlapped with gene promoter regions
# loop regions a
bedtools intersect -a ${od}/8cell_significant_loops.bed -b ${pro1k} -F 0.5 -wa > ${od}/loops_a_overlapped_pro1k.bed
# loop regions b
cat ${od}/8cell_significant_loops.bed | awk 'BEGIN{FS=OFS="\t"}{print $4,$5,$6,$1,$2,$3,$7,$8}' \
    | bedtools intersect -a stdin -b ${pro1k} -F 0.5 -wa \
    | awk 'BEGIN{FS=OFS="\t"}{print $4,$5,$6,$1,$2,$3,$7,$8}' > ${od}/loops_b_overlapped_pro1k.bed
# merge and remove duplicated lines
loops_bed=${od}/8cell_significant_loops_overlapped_with_pro1k.bed
cat ${od}/loops_a_overlapped_pro1k.bed ${od}/loops_b_overlapped_pro1k.bed | sort -u > ${loops_bed}
rm ${od}/loops_a* ${od}/loops_b*

# - ChIPseq&ATACseq + HiCsubset distal H3K27ac and ATAC peaks overlapped with significant loops regions
# merge loops regions into one bed file
loops_region=${od}/8cell_significant_loops_overlapped_with_pro1k_regions.bed
awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3"\n"$4,$5,$6}' ${loops_bed} | sort -u > ${loops_region}
# subset distal H3K27ac and ATAC peaks overlapped with significant loops regions
c8_enhancer=${od}/8cell_enhancer.bed
cat ${od}/8cell_distal_h3k27ac_pks.bed ${od}/8cell_distal_atac_pks.bed \
    | sort -k1,1 -k2,2n \
    | bedtools intersect -a stdin -b ${loops_region} -f 0.5 -wa > ${c8_enhancer}

# - ChIPseq&ATACseq: define reprogrammed enhancer and primed enhancer
c8_enhancer=${od}/8cell_enhancer.bed
c8_RE_k93=${od}/8cell_reprogrammed_enhancer_H3K9me3.bed
c8_RE_k273=${od}/8cell_reprogrammed_enhancer_H3K27me3.bed
c8_PE=${od}/8cell_primed_enhancer.bed
# - H3K9me3-reprogrammed
bedtools intersect -a ${c8_enhancer} -b ${c4_93_pks} -f 0.5 \
   | awk '{if(($3-$2)>=200) print $0}' | sort -u > ${c8_RE_k93}
# - H3K27me3-reprogrammed
bedtools intersect -a ${c8_enhancer} -b ${c4_273_pks} -f 0.5 -wa \
   | bedtools intersect -a stdin -b ${c2_273_pks} -f 0.5 \
   | awk '{if(($3-$2)>=200) print $0}' | sort -u > ${c8_RE_k273}
# - primed
bedtools intersect -v -a ${c8_enhancer} -b ${c4_93_pks} \
   | bedtools intersect -a stdin -b ${c4_atac_pks} -f 0.5 \
   | awk '{if(($3-$2)>=200) print $0}' | sort -u > ${c8_PE}
# SVA-derived enhancers
c8_enhancer=${od}/8cell_enhancer.bed
c8_RE_k93=${od}/8cell_reprogrammed_enhancer_H3K9me3.bed
c8_RE_k273=${od}/8cell_reprogrammed_enhancer_H3K27me3.bed
c8_PE=${od}/8cell_primed_enhancer.bed
sva_c8_enhancer=${od}/8cell_enhancer_derived_by_SVA.bed
sva_c8_RE_k93=${od}/8cell_reprogrammed_enhancer_H3K9me3_derived_by_SVA.bed
sva_c8_RE_k273=${od}/8cell_reprogrammed_enhancer_H3K27me3_derived_by_SVA.bed
sva_c8_RE=${od}/8cell_reprogrammed_enhancer_derived_by_SVA.bed
sva_c8_PE=${od}/8cell_primed_enhancer_derived_by_SVA.bed
repeat_anno=/home/yhw/genome/ucsc/homo_sapiens/hg38/repeat/GRCh38_RepeatMasker_repeat_region.bed
cat ${repeat_anno} | grep "SVA" | bedtools intersect -a ${c8_enhancer} -b stdin -wa -wb | awk '$1 !~ "_"'> ${sva_c8_enhancer}
cat ${repeat_anno} | grep "SVA" | bedtools intersect -a ${c8_RE_k93} -b stdin -wa -wb > ${sva_c8_RE_k93}
cat ${repeat_anno} | grep "SVA" | bedtools intersect -a ${c8_RE_k273} -b stdin -wa -wb > ${sva_c8_RE_k273}
cat ${repeat_anno} | grep "SVA" | bedtools intersect -a ${c8_PE} -b stdin -wa -wb > ${sva_c8_PE}
cat ${sva_c8_RE_k273} ${sva_c8_RE_k93} > ${sva_c8_RE}
sva_no_enhancer=${od}/SVAs_not_in_8cell_enhancer.bed
cut -f 7 ${sva_c8_enhancer} | grep -w -v -f - ${repeat_anno} | grep "SVA" > ${sva_no_enhancer}



### >>> Plotting signal on PE and RE (YHW)
# - region
c8_RE_k93=${od}/8cell_reprogrammed_enhancer_H3K9me3.bed
c8_RE_k273=${od}/8cell_reprogrammed_enhancer_H3K27me3.bed
c8_PE=${od}/8cell_primed_enhancer.bed
# - bigwig files
# H3K27ac
c8_27ac_1_bw=/home/yhw/bioinfo/project-temp/wjc/gse124718/results/deeptools/coverage/SRR9131735_1_val_1_SRR9131735_2_val_2_psort_dedup.bw
c8_27ac_2_bw=/home/yhw/bioinfo/project-temp/wjc/gse124718/results/deeptools/coverage/SRR9131736_1_val_1_SRR9131736_2_val_2_psort_dedup.bw
c8_27ac_bw=/home/yhw/bioinfo/project-temp/wjc/gse124718/results/deeptools/coverage/merge/8cell_merge.bw
# H3K9me3
c4_93_bw=/home/cmq/bioinfo/project-cmq/embryo_93/results/deeptools/bigwig/merged_bam_ratio/h4C93-Ip-1_vs_h4C93-Input-1.bw
c8_93_bw=/home/cmq/bioinfo/project-cmq/embryo_93/results/deeptools/bigwig/merged_bam_ratio/h8C93-Ip-3_vs_h8C93-Input-2.bw
# H3K27me3
c2_273_bw=/home/cmq/bioinfo/project-cmq/gse124718/results/deeptools/bigwig/merge_replicate_multi/2cell_multi_replicate_merged.bw
c4_273_bw=/home/cmq/bioinfo/project-cmq/gse124718/results/deeptools/bigwig/merge_replicate_multi/4cell_multi_replicate_merged.bw
c8_273_bw=/home/cmq/bioinfo/project-cmq/gse124718/results/deeptools/bigwig/merge_replicate_multi/8cell_multi_replicate_merged.bw
# ATAC
atac_dir=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/deeptools/bigwig/dedup_merge_nfr_bam
zy_atac_1_bw=${atac_dir}/SRR6501974_1_val_1_SRR6501974_2_val_2_psorted_dedup.bw
c2_atac_2_bw=${atac_dir}/SRR5837317_1_val_1_SRR5837317_2_val_2_psorted_dedup_SRR5837318_1_val_1_SRR5837318_2_val_2_psorted_dedup_merged.bw
c4_atac_2_bw=${atac_dir}/SRR5837326_1_val_1_SRR5837326_2_val_2_psorted_dedup_SRR5837327_1_val_1_SRR5837327_2_val_2_psorted_dedup_SRR5837328_1_val_1_SRR5837328_2_val_2_psorted_dedup_merged.bw
atac_dir=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/deeptools/bigwig/dedup_merge_bam
zy_atac_2_bw=${atac_dir}/SRR6501974_1_val_1_SRR6501974_2_val_2_psorted_dedup.bw
c2_atac_1_bw=${atac_dir}/SRR5837315_1_val_1_SRR5837315_2_val_2_psorted_dedup_SRR5837316_1_val_1_SRR5837316_2_val_2_psorted_dedup_merged.bw
c4_atac_1_bw=${atac_dir}/SRR5837324_1_val_1_SRR5837324_2_val_2_psorted_dedup_SRR5837325_1_val_1_SRR5837325_2_val_2_psorted_dedup_merged.bw
c8_atac_1_bw=${atac_dir}/SRR5837329_1_val_1_SRR5837329_2_val_2_psorted_dedup_SRR5837330_1_val_1_SRR5837330_2_val_2_psorted_dedup_merged.bw
c8_atac_2_bw=${atac_dir}/SRR5837331_1_val_1_SRR5837331_2_val_2_psorted_dedup_SRR5837332_1_val_1_SRR5837332_2_val_2_psorted_dedup_merged.bw
# DNase
c2_dna_1_bw=/home/data/publicdata/PRJCA000484/analysis/results/deeptools/coverage/CRR019752_f1_val_1_CRR019752_r2_val_2_psort_dedup.bw
c2_dna_2_bw=/home/data/publicdata/PRJCA000484/analysis/results/deeptools/coverage/CRR019753_f1_val_1_CRR019753_r2_val_2_psort_dedup.bw
c4_dna_1_bw=/home/data/publicdata/PRJCA000484/analysis/results/deeptools/coverage/CRR019754_f1_val_1_CRR019754_r2_val_2_psort_dedup.bw
c4_dna_2_bw=/home/data/publicdata/PRJCA000484/analysis/results/deeptools/coverage/CRR019755_f1_val_1_CRR019755_r2_val_2_psort_dedup.bw
c8_dna_1_bw=/home/data/publicdata/PRJCA000484/analysis/results/deeptools/coverage/CRR019756_f1_val_1_CRR019756_r2_val_2_psort_dedup.bw
c8_dna_2_bw=/home/data/publicdata/PRJCA000484/analysis/results/deeptools/coverage/CRR019757_f1_val_1_CRR019757_r2_val_2_psort_dedup.bw
mkdir -p /home/cmq/bioinfo/project-cmq/embryo_93/results/deeptools/matrix/revised_8c_enhancer
mkdir -p /home/cmq/bioinfo/project-cmq/embryo_93/results/deeptools/heatmap/revised_8c_enhancer
computeMatrix reference-point -S ${c8_27ac_bw} ${c8_27ac_1_bw} ${c8_27ac_2_bw} ${c2_273_bw} ${c4_273_bw} ${c8_273_bw} ${c4_93_bw} ${c8_93_bw} \
                                 ${zy_atac_1_bw} ${zy_atac_2_bw} ${c2_atac_1_bw} ${c2_atac_2_bw} ${c4_atac_1_bw} ${c4_atac_2_bw} ${c8_atac_1_bw} ${c8_atac_2_bw} \
                                 ${c2_dna_1_bw} ${c2_dna_2_bw} ${c4_dna_1_bw} ${c4_dna_2_bw} ${c8_dna_1_bw} ${c8_dna_2_bw} \
                              -R ${c8_RE_k93} ${c8_RE_k273} ${c8_PE} --averageTypeBins mean --referencePoint center --beforeRegionStartLength 3000 --afterRegionStartLength 3000 \
                              --numberOfProcessors 30 -o ../results/deeptools/matrix/revised_8c_enhancer/enhancer_gzmatrix.gz
# - H3K27ac signal: Greens
for max in `seq 10 5 50`
do
    plotHeatmap -m ../results/deeptools/matrix/revised_8c_enhancer/enhancer_gzmatrix.gz \
                -o ../results/deeptools/heatmap/revised_8c_enhancer/different_signal_on_8cell_enhancer_H3K27ac_heatmap_rpkm_0-${max}.pdf \
                --averageTypeSummaryPlot mean  --heatmapHeight 8 --heatmapWidth 2 --whatToShow "plot, heatmap and colorbar" \
                --plotType lines --colorMap Greens --zMin 0 --zMax ${max} --yMin 0 --yMax 25 \
                --legendLocation upper-left --interpolationMethod auto --missingDataColor 1 \
                --samplesLabel "8c 27ac" "8c 27ac 1" "8c 27ac 2" "2c 273" "4c 273" "8c 273" "4c 93" "8c 93" "zygote atac1" "zygote atac2" "2c atac1" "2c atac2" "4c atac1" "4c atac2" "8c atac1" "8c atac2" "2c dna1" "2c dna2" "4c dna1" "4c dna2" "8c dna1" "8c dna2" \
                --regionsLabel "RE.93" "RE.273" "PE"
done
# - H3K27me3 signal: Oranges
for max in `seq 10 5 50`
do
    plotHeatmap -m ../results/deeptools/matrix/revised_8c_enhancer/enhancer_gzmatrix.gz \
                -o ../results/deeptools/heatmap/revised_8c_enhancer/different_signal_on_8cell_enhancer_H3K27me3_heatmap_rpkm_0-${max}.pdf \
                --averageTypeSummaryPlot mean  --heatmapHeight 8 --heatmapWidth 2 --whatToShow "plot, heatmap and colorbar" \
                --plotType lines --colorMap Oranges --zMin 0 --zMax ${max} --yMin 0 --yMax 50 \
                --legendLocation upper-left --interpolationMethod auto --missingDataColor 1 \
                --samplesLabel "8c 27ac" "8c 27ac 1" "8c 27ac 2" "2c 273" "4c 273" "8c 273" "4c 93" "8c 93" "zygote atac1" "zygote atac2" "2c atac1" "2c atac2" "4c atac1" "4c atac2" "8c atac1" "8c atac2" "2c dna1" "2c dna2" "4c dna1" "4c dna2" "8c dna1" "8c dna2" \
                --regionsLabel "RE.93" "RE.273" "PE"
done
# - ATAC signal: Blues
for max in `seq 10 5 50`
do
    plotHeatmap -m ../results/deeptools/matrix/revised_8c_enhancer/enhancer_gzmatrix.gz \
                -o ../results/deeptools/heatmap/revised_8c_enhancer/different_signal_on_8cell_enhancer_ATAC_heatmap_rpkm_0-${max}.pdf \
                --averageTypeSummaryPlot mean  --heatmapHeight 8 --heatmapWidth 2 --whatToShow "plot, heatmap and colorbar" \
                --plotType lines --colorMap Blues --zMin 0 --zMax ${max} --yMin 0 --yMax 80 \
                --legendLocation upper-left --interpolationMethod auto --missingDataColor 1 \
                --samplesLabel "8c 27ac" "8c 27ac 1" "8c 27ac 2" "2c 273" "4c 273" "8c 273" "4c 93" "8c 93" "zygote atac1" "zygote atac2" "2c atac1" "2c atac2" "4c atac1" "4c atac2" "8c atac1" "8c atac2" "2c dna1" "2c dna2" "4c dna1" "4c dna2" "8c dna1" "8c dna2" \
                --regionsLabel "RE.93" "RE.273" "PE"
done
# - H3K9me3 signal: cyan, white, red
plotHeatmap -m ../results/deeptools/matrix/revised_8c_enhancer/enhancer_gzmatrix.gz \
            -o ../results/deeptools/heatmap/revised_8c_enhancer/different_signal_on_8cell_enhancer_H3K9me3_heatmap_lfc.pdf \
            --averageTypeSummaryPlot mean  --heatmapHeight 8 --heatmapWidth 2 --whatToShow "plot, heatmap and colorbar" \
            --plotType lines --colorList cyan,white,red --zMin -1 --zMax 2 --yMin -1 --yMax 2 \
            --legendLocation upper-left --interpolationMethod auto --missingDataColor 1 \
            --samplesLabel "8c 27ac" "8c 27ac 1" "8c 27ac 2" "2c 273" "4c 273" "8c 273" "4c 93" "8c 93" "zygote atac1" "zygote atac2" "2c atac1" "2c atac2" "4c atac1" "4c atac2" "8c atac1" "8c atac2" "2c dna1" "2c dna2" "4c dna1" "4c dna2" "8c dna1" "8c dna2" \
            --regionsLabel "RE.93" "RE.273" "PE"



### >>> Footprint analysis (YHW)
wd=/home/cmq/bioinfo/project-cmq/embryo_93
bam_dir=/home/cmq/bioinfo/project-cmq/embryo_93/results/samtools/merge_atac
[ ! -e "${bam_dir}" ] && mkdir -p ${bam_dir}

# - load peaks
c8_RE_k93=${od}/8cell_reprogrammed_enhancer_H3K9me3.bed
c8_RE_k273=${od}/8cell_reprogrammed_enhancer_H3K27me3.bed
c8_PE=${od}/8cell_primed_enhancer.bed
sva_c8_enhancer=${od}/8cell_enhancer_derived_by_SVA.bed
sva_c8_RE_k93=${od}/8cell_reprogrammed_enhancer_H3K9me3_derived_by_SVA.bed
sva_c8_RE_k273=${od}/8cell_reprogrammed_enhancer_H3K27me3_derived_by_SVA.bed
sva_c8_PE=${od}/8cell_primed_enhancer_derived_by_SVA.bed
# - load files for TOBIAS
jaspar=/home/yhw/bioinfo/project-mine/Embryo.93/chipseq/hs_k9me3/metadata/motif.jaspar
blacklist=/home/yhw/software/blender/hg38.blacklist.bed
genome_fa=/home/yhw/genome/ucsc/homo_sapiens/hg38/hg38.fa

# - load bam files
atac_dir=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/samtools/merge/dedup_bam
zy_atac_bam=${atac_dir}/SRR6501974_1_val_1_SRR6501974_2_val_2_psorted_dedup.bam
c2_atac_bam=${atac_dir}/SRR5837315_1_val_1_SRR5837315_2_val_2_psorted_dedup_SRR5837316_1_val_1_SRR5837316_2_val_2_psorted_dedup_merged.bam
c4_atac_bam=${atac_dir}/SRR5837324_1_val_1_SRR5837324_2_val_2_psorted_dedup_SRR5837325_1_val_1_SRR5837325_2_val_2_psorted_dedup_merged.bam
c8_atac_bam=${atac_dir}/SRR5837329_1_val_1_SRR5837329_2_val_2_psorted_dedup_SRR5837330_1_val_1_SRR5837330_2_val_2_psorted_dedup_merged.bam
cp ${zy_atac_bam} ${bam_dir}/zygote.bam
cp ${c2_atac_bam} ${bam_dir}/2cell.bam
cp ${c4_atac_bam} ${bam_dir}/4cell.bam
cp ${c8_atac_bam} ${bam_dir}/8cell.bam
ls ${bam_dir} | grep "bam$" > ${wd}/metadata/footprint_bam_list.txt

# - run TOBIAS
FPtobias(){
   indir=$1; file=$2; outdir=$3; logdir=$4; bed=$5; seq=$6; motif=$7; bl=$8
   # - indir: input dir includes the bam files;
   # - file: txt file records the names of bam files;
   # - outdir: output directory;
   # - logdir: log file directory;
   # - seq: genome sequence fa file;
   # - blacklist: bed file record the black list regions
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   check_index=$(ls ${indir}/* | grep "bai$")
   if [ "${check_index}" -eq 0 ]; then
      find ${bam_dir} -name "*.bam" | parallel -P 5 --gnu "samtools index -@ 12 -b {}"
   fi
   bedname=$(basename -s ".bed" ${bed})
   outdir=${outdir}/${bedname}; logdir=${logdir}/${bedname}
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   # - ATACorrect
   find ${bam_dir} -name "*.bam" | xargs basename -s ".bam" \
      | parallel -P 3 --gnu "TOBIAS ATACorrect --bam ${indir}/{}.bam --genome ${seq} --peaks ${bed} --blacklist ${bl} --outdir ${outdir} --cores 10 --prefix {} > ${logdir}/{}.log 2>&1"
   # - ScoreBigwig
   find ${outdir} -name "*_corrected.bw" | xargs basename -s ".bw" \
      | parallel -P 3 --gnu "TOBIAS ScoreBigwig --signal ${outdir}/{}.bw --regions ${bed} --output ${outdir}/{}_footprint.bw --cores 10 > ${logdir}/{}.log 2>&1"
   # - BINDetect
   fbw_list=$(ls ${outdir}/*_footprint.bw | tr "\n" " ")
   bw_name=$(ls ${indir}/*bam | xargs basename -s ".bam" | tr "\n" " ")
   TOBIAS BINDetect --motifs ${motif} --signals ${fbw_list} --genome ${seq} --peaks ${bed} \
                    --outdir ${outdir}/BIND_output --cond_names ${bw_name} --cores 20 > ${logdir}/BIND.log 2>&1
   # - PlotAggregate
   [ ! -d "${outdir}/PlotAggregate_output" ] && mkdir -p ${outdir}/PlotAggregate_output
   find ${outdir}/BIND_output/ -name "*" -type d -maxdepth 1 | awk 'NR>1' | while read dir
   do
       bed_list=$(ls ${dir}/beds/*bed | grep "bed$" | tr "\n" " ")
       echo ${bed_list}
       cbw_list=$(ls ${outdir}/*_corrected.bw | tr "\n" " ")
       TOBIAS PlotAggregate --TFBS ${bed_list} --signals ${cbw_list} \
                            --output ${outdir}/PlotAggregate_output/${dir##*/}_footprint_comparison.pdf \
                            --share_y both --plot_boundaries --signal-on-x --smooth 8 --flank 50
   done
}
# - PEs
FPtobias ${bam_dir} ${wd}/metadata/footprint_bam_list.txt ${wd}/results/tobias ${wd}/logs/tobias \
         ${c8_PE} ${genome_fa} ${jaspar} ${blacklist}
# - REs H3K9me3
FPtobias ${bam_dir} ${wd}/metadata/footprint_bam_list.txt ${wd}/results/tobias ${wd}/logs/tobias \
         ${c8_RE_k93} ${genome_fa} ${jaspar} ${blacklist}
# - REs H3K27me3
FPtobias ${bam_dir} ${wd}/metadata/footprint_bam_list.txt ${wd}/results/tobias ${wd}/logs/tobias \
         ${c8_RE_k273} ${genome_fa} ${jaspar} ${blacklist}
# - SVA PEs
FPtobias ${bam_dir} ${wd}/metadata/footprint_bam_list.txt ${wd}/results/tobias ${wd}/logs/tobias \
         ${sva_c8_PE} ${genome_fa} ${jaspar} ${blacklist}
# - SVA REs H3K9me3
FPtobias ${bam_dir} ${wd}/metadata/footprint_bam_list.txt ${wd}/results/tobias ${wd}/logs/tobias \
         ${sva_c8_RE_k93} ${genome_fa} ${jaspar} ${blacklist}
# - SVA REs H3K27me3
FPtobias ${bam_dir} ${wd}/metadata/footprint_bam_list.txt ${wd}/results/tobias ${wd}/logs/tobias \
         ${sva_c8_RE_k273} ${genome_fa} ${jaspar} ${blacklist}
# - Find proteins interacted with K9-factor proteins
# K9 factor list: SETDB1, SETDB2, SUV39H1, SUV39H2, UHRF1, CBX1, CBX5, TRIM28
hs_biogrid=/home/yhw/document/biogrid/4.4.204/BIOGRID-ORGANISM-Homo_sapiens-4.4.204.tab3.txt
hs_k9_factor=../metadata/K9_factor_list.txt
hs_k9_factor_binds=../metadata/K9_factor_list_binding_proteins.txt
cat ${hs_biogrid} | sed 's,^#,,g' | awk 'BEGIN{FS=OFS="\t"}{print $8,$9,$12,$13,$15}' | grep -w -f ${hs_k9_factor} - > ${hs_k9_factor_binds}
hs_k9_eraser=../metadata/K9_erasers_list.txt
hs_k9_eraser_binds=../metadata/K9_erasers_list_binding_proteins.txt
cat ${hs_biogrid} | sed 's,^#,,g' | awk 'BEGIN{FS=OFS="\t"}{print $8,$9,$12,$13,$15}' | grep -w -f ${hs_k9_eraser} - > ${hs_k9_eraser_binds}



### >>> Find the genes interacted by self-defined enhancers: RE.K9, RE.K273, PE (YHW)
od=/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer
[ ! -e "${od}" ] && mkdir -p ${od}

# - load gene promoter regions
hs_pro1k=/home/yhw/genome/ensembl/release97/homo_sapiens/dna_anno/regions/Homo_sapiens.GRCh38.97_pro1k_ucsc.bed

# - load enhancers
c8_RE_k93=${od}/8cell_reprogrammed_enhancer_H3K9me3.bed
c8_RE_k273=${od}/8cell_reprogrammed_enhancer_H3K27me3.bed
c8_PE=${od}/8cell_primed_enhancer.bed

# - load interaction loops
c2_loops=/home/yhw/bioinfo/project-mine/Embryo.93/hic/hs_embryo/results/fithic/2cell_rep1_loops/FitHiC.spline_pass1.res40000.significances.txt.gz
c8_loops=/home/yhw/bioinfo/project-mine/Embryo.93/hic/hs_embryo/results/fithic/8cell_rep2_loops/FitHiC.spline_pass1.res40000.significances.txt.gz

# - 2cell loops with genes
for enhancer in c8_RE_k93 c8_RE_k273 c8_PE
do
    zcat ${c2_loops} \
       | awk '{OFS=FS="\t"}{if($7<=0.005&&$5>=30) print $1,$2-20000,$2+20000,$3,$4-20000,$4+20000,$5,$6,$7,$8,$9}' \
       | awk '{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' \
       | bedtools intersect -a stdin -b ${!enhancer} -wa -F 1 | cut -f 4-6 \
       | bedtools intersect -a stdin -b ${hs_pro1k} -wa -wb -F 1 \
       | cut -f 4-8 > ${od}/genes_interacted_with_${enhancer}_enhancers_in_2cell.txt
    zcat ${c2_loops} \
       | awk '{OFS=FS="\t"}{if($7<=0.05&&$5>=30) print $1,$2-20000,$2+20000,$3,$4-20000,$4+20000,$5,$6,$7,$8,$9}' \
       | awk '{FS=OFS="\t"}{print $4,$5,$6,$1,$2,$3,$7,$8,$9}' \
       | bedtools intersect -a stdin -b ${!enhancer} -wa -F 1 | cut -f 4-6 \
       | bedtools intersect -a stdin -b ${hs_pro1k} -wa -wb -F 1 \
       | cut -f 4-8 >> ${od}/genes_interacted_with_${enhancer}_enhancers_in_2cell.txt
    sort -u ${od}/genes_interacted_with_${enhancer}_enhancers_in_2cell.txt | sponge ${od}/genes_interacted_with_${enhancer}_enhancers_in_2cell.txt
done
# - 8cell loops with genes
for enhancer in c8_RE_k93 c8_RE_k273 c8_PE
do
    zcat ${c8_loops} \
       | awk '{OFS=FS="\t"}{if($7<=0.05&&$5>=30) print $1,$2-20000,$2+20000,$3,$4-20000,$4+20000,$5,$6,$7,$8,$9}' \
       | awk '{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' \
       | bedtools intersect -a stdin -b ${!enhancer} -wa -F 1 | cut -f 4-6 \
       | bedtools intersect -a stdin -b ${hs_pro1k} -wa -wb -F 1 \
       | cut -f 4-8 > ${od}/genes_interacted_with_${enhancer}_enhancers_in_8cell.txt
    zcat ${c8_loops} \
       | awk '{OFS=FS="\t"}{if($7<=0.05&&$5>=30) print $1,$2-20000,$2+20000,$3,$4-20000,$4+20000,$5,$6,$7,$8,$9}' \
       | awk '{FS=OFS="\t"}{print $4,$5,$6,$1,$2,$3,$7,$8,$9}' \
       | bedtools intersect -a stdin -b ${!enhancer} -wa -F 1 | cut -f 4-6 \
       | bedtools intersect -a stdin -b ${hs_pro1k} -wa -wb -F 1 \
       | cut -f 4-8 >> ${od}/genes_interacted_with_${enhancer}_enhancers_in_8cell.txt
    sort -u ${od}/genes_interacted_with_${enhancer}_enhancers_in_8cell.txt | sponge ${od}/genes_interacted_with_${enhancer}_enhancers_in_8cell.txt
done



### >>> Find the genes interacted by self-defined enhancers: RE.K9, RE.K273, PE (YHW)
od=/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer
[ ! -e "${od}" ] && mkdir -p ${od}

# - load gene promoter regions
hs_pro1k=/home/yhw/genome/ensembl/release97/homo_sapiens/dna_anno/regions/Homo_sapiens.GRCh38.97_pro1k_ucsc.bed

# - load enhancers
c8_RE_k93=${od}/8cell_reprogrammed_enhancer_H3K9me3.bed
c8_RE_k273=${od}/8cell_reprogrammed_enhancer_H3K27me3.bed
c8_PE=${od}/8cell_primed_enhancer.bed
sva_c8_enhancer=${od}/8cell_enhancer_derived_by_SVA.bed
sva_c8_RE_k93=${od}/8cell_reprogrammed_enhancer_H3K9me3_derived_by_SVA.bed
sva_c8_RE_k273=${od}/8cell_reprogrammed_enhancer_H3K27me3_derived_by_SVA.bed
sva_c8_PE=${od}/8cell_primed_enhancer_derived_by_SVA.bed

# - load interaction loops
c2_loops=/home/yhw/bioinfo/project-mine/Embryo.93/hic/hs_embryo/results/fithic/2cell_rep1_loops/FitHiC.spline_pass1.res40000.significances.txt.gz
c8_loops=/home/yhw/bioinfo/project-mine/Embryo.93/hic/hs_embryo/results/fithic/8cell_rep2_loops/FitHiC.spline_pass1.res40000.significances.txt.gz

# - 2cell loops with genes
for enhancer in c8_RE_k93 c8_RE_k273 c8_PE
do
    zcat ${c2_loops} \
       | awk '{OFS=FS="\t"}{if($7<=0.005&&$5>=30) print $1,$2-20000,$2+20000,$3,$4-20000,$4+20000,$5,$6,$7,$8,$9}' \
       | awk '{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' \
       | bedtools intersect -a stdin -b ${!enhancer} -wa -wb -F 1 \
       | awk '{FS=OFS="\t"}{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12}' \
       | bedtools intersect -a stdin -b ${hs_pro1k} -wa -wb -F 1 \
       | awk '{FS=OFS="\t"}{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' > ${od}/genes_interacted_with_${enhancer}_enhancers_in_2cell_full.txt
    zcat ${c2_loops} \
       | awk '{OFS=FS="\t"}{if($7<=0.005&&$5>=30) print $1,$2-20000,$2+20000,$3,$4-20000,$4+20000,$5,$6,$7,$8,$9}' \
       | awk '{FS=OFS="\t"}{print $4,$5,$6,$1,$2,$3,$7,$8,$9}' \
       | bedtools intersect -a stdin -b ${!enhancer} -wa -wb -F 1 \
       | awk '{FS=OFS="\t"}{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12}' \
       | bedtools intersect -a stdin -b ${hs_pro1k} -wa -wb -F 1 \
       | awk '{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' >> ${od}/genes_interacted_with_${enhancer}_enhancers_in_2cell_full.txt
    sort -u ${od}/genes_interacted_with_${enhancer}_enhancers_in_2cell_full.txt | sponge ${od}/genes_interacted_with_${enhancer}_enhancers_in_2cell_full.txt
done
for enhancer in sva_c8_RE_k93 sva_c8_RE_k273 sva_c8_PE
do
    zcat ${c2_loops} \
       | awk '{OFS=FS="\t"}{if($6<=0.05&&$5>=20) print $1,$2-20000,$2+20000,$3,$4-20000,$4+20000,$5,$6,$7,$8,$9}' \
       | awk '{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' \
       | bedtools intersect -a stdin -b ${!enhancer} -wa -wb -F 1 \
       | awk '{FS=OFS="\t"}{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12}' \
       | bedtools intersect -a stdin -b ${hs_pro1k} -wa -wb -F 1 \
       | awk '{FS=OFS="\t"}{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' > ${od}/genes_interacted_with_${enhancer}_enhancers_in_2cell_full.txt
    zcat ${c2_loops} \
       | awk '{OFS=FS="\t"}{if($6<=0.05&&$5>=20) print $1,$2-20000,$2+20000,$3,$4-20000,$4+20000,$5,$6,$7,$8,$9}' \
       | awk '{FS=OFS="\t"}{print $4,$5,$6,$1,$2,$3,$7,$8,$9}' \
       | bedtools intersect -a stdin -b ${!enhancer} -wa -wb -F 1 \
       | awk '{FS=OFS="\t"}{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12}' \
       | bedtools intersect -a stdin -b ${hs_pro1k} -wa -wb -F 1 \
       | awk '{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' >> ${od}/genes_interacted_with_${enhancer}_enhancers_in_2cell_full.txt
    sort -u ${od}/genes_interacted_with_${enhancer}_enhancers_in_2cell_full.txt | sponge ${od}/genes_interacted_with_${enhancer}_enhancers_in_2cell_full.txt
done

# - 8cell loops with genes
for enhancer in c8_RE_k93 c8_RE_k273 c8_PE
do
    zcat ${c8_loops} \
       | awk '{OFS=FS="\t"}{if($7<=0.05&&$5>=30) print $1,$2-20000,$2+20000,$3,$4-20000,$4+20000,$5,$6,$7,$8,$9}' \
       | awk '{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' \
       | bedtools intersect -a stdin -b ${!enhancer} -wa -wb -F 1 \
       | awk '{FS=OFS="\t"}{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12}' \
       | bedtools intersect -a stdin -b ${hs_pro1k} -wa -wb -F 1 \
       | awk '{FS=OFS="\t"}{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' > ${od}/genes_interacted_with_${enhancer}_enhancers_in_8cell_full.txt
    zcat ${c8_loops} \
       | awk '{OFS=FS="\t"}{if($7<=0.05&&$5>=30) print $1,$2-20000,$2+20000,$3,$4-20000,$4+20000,$5,$6,$7,$8,$9}' \
       | awk '{FS=OFS="\t"}{print $4,$5,$6,$1,$2,$3,$7,$8,$9}' \
       | bedtools intersect -a stdin -b ${!enhancer} -wa -wb -F 1 \
       | awk '{FS=OFS="\t"}{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12}' \
       | bedtools intersect -a stdin -b ${hs_pro1k} -wa -wb -F 1 \
       | awk '{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' >> ${od}/genes_interacted_with_${enhancer}_enhancers_in_8cell_full.txt
    sort -u ${od}/genes_interacted_with_${enhancer}_enhancers_in_8cell_full.txt | sponge ${od}/genes_interacted_with_${enhancer}_enhancers_in_8cell_full.txt
done
for enhancer in sva_c8_RE_k93 sva_c8_RE_k273 sva_c8_PE
do
    zcat ${c8_loops} \
       | awk '{OFS=FS="\t"}{if($6<=0.05&&$5>=20) print $1,$2-20000,$2+20000,$3,$4-20000,$4+20000,$5,$6,$7,$8,$9}' \
       | awk '{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' \
       | bedtools intersect -a stdin -b ${!enhancer} -wa -wb -F 1 \
       | awk '{FS=OFS="\t"}{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12}' \
       | bedtools intersect -a stdin -b ${hs_pro1k} -wa -wb -F 1 \
       | awk '{FS=OFS="\t"}{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' > ${od}/genes_interacted_with_${enhancer}_enhancers_in_8cell_full.txt
    zcat ${c8_loops} \
       | awk '{OFS=FS="\t"}{if($6<=0.05&&$5>=20) print $1,$2-20000,$2+20000,$3,$4-20000,$4+20000,$5,$6,$7,$8,$9}' \
       | awk '{FS=OFS="\t"}{print $4,$5,$6,$1,$2,$3,$7,$8,$9}' \
       | bedtools intersect -a stdin -b ${!enhancer} -wa -wb -F 1 \
       | awk '{FS=OFS="\t"}{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12}' \
       | bedtools intersect -a stdin -b ${hs_pro1k} -wa -wb -F 1 \
       | awk '{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' >> ${od}/genes_interacted_with_${enhancer}_enhancers_in_8cell_full.txt
    sort -u ${od}/genes_interacted_with_${enhancer}_enhancers_in_8cell_full.txt | sponge ${od}/genes_interacted_with_${enhancer}_enhancers_in_8cell_full.txt
done

# - Annotate loops with repeat
od=/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer
repeat_anno=/home/yhw/genome/ucsc/homo_sapiens/hg38/repeat/GRCh38_RepeatMasker_repeat_region.bed
ls ../results/bedtools/revised_8c_enhancer/*8cell_full.txt | grep "sva" | while read file
do
    prefix=$(basename -s ".txt" ${file})
    echo ${prefix}
    awk '{FS=OFS="\t"}{print $10,$11,$12,$13,$14,$15,$16,$17}' ${file} \
       | bedtools intersect -a stdin -b ${repeat_anno} -f 0.25 -wa -wb | grep "SVA" > ${od}/${prefix}_anno_by_repeat.txt
    awk '{FS=OFS="\t"}{if(($5-$2)>200000) print $10,$11,$12,$13,$14,$15,$16,$17,$1,$2,$3,$4,$5,$6,$7,$8,$9}' ${file} \
       | bedtools intersect -a stdin -b ${repeat_anno} -f 0.25 -wa -wb | grep "SVA" > ${od}/${prefix}_anno_by_repeat_distal_longer_4bins.txt
done



### >>> Interaction aggregate plots (YHW)
# - output dir
od=/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer
# - process file for GREAT
c8_RE_k93=${od}/8cell_reprogrammed_enhancer_H3K9me3.bed
c8_RE_k273=${od}/8cell_reprogrammed_enhancer_H3K27me3.bed
c8_PE=${od}/8cell_primed_enhancer.bed
awk '{print $0"\t"$1":"$2"-"$3}' ${c8_RE_k93} > ${od}/8cell_reprogrammed_enhancer_H3K9me3_for_GREAT.bed
awk '{print $0"\t"$1":"$2"-"$3}' ${c8_RE_k273} > ${od}/8cell_reprogrammed_enhancer_H3K27me3_for_GREAT.bed
awk '{print $0"\t"$1":"$2"-"$3}' ${c8_PE} > ${od}/8cell_primed_enhancer_for_GREAT.bed

# - process gene-region txt file
hs_pro3k=/home/yhw/genome/ensembl/release97/homo_sapiens/dna_anno/regions/Homo_sapiens.GRCh38.97_pro3k_ucsc.bed
sed 's,:,\t,g' ${hs_pro3k} | sort -k5,5 > ../metadata/Homo_sapiens.GRCh38.97_pro3k_ucsc.bed
less ../results/bedtools/revised_8c_enhancer/8cell_reprogrammed_enhancer_H3K9me3_for_GREAT_region_gene_table.txt \
   | sed -e 's, (,\t,g' | sed 's,)\,,\t,g' | sed 's,),,g' | sed 's,\t ,\t,g' \
   | awk '{print $1"\t"$2"\n"$1"\t"$4"\n"$1"\t"$6"\n"$1"\t"$8}' | awk '{if($2!="") print $0}' \
   | grep -v "^#" | sort -u | sort -k2,2 \
   | join -1 2 -2 5 - ../metadata/Homo_sapiens.GRCh38.97_pro3k_ucsc.bed \
   | sed 's, ,\t,g' | cut -f 2-5 \
   | sed -e 's,:,\t,g' -e 's,-,\t,g' > ../results/bedtools/revised_8c_enhancer/8cell_reprogrammed_enhancer_H3K9me3_for_GREAT_region_gene_table.bed
c8_RE_k93_great_bedpe=../results/bedtools/revised_8c_enhancer/8cell_reprogrammed_enhancer_H3K9me3_for_GREAT_region_gene_table.bedpe
cut -f 1-6 ../results/bedtools/revised_8c_enhancer/8cell_reprogrammed_enhancer_H3K9me3_for_GREAT_region_gene_table.bed | sort -k1,1 -k2,2n > $c8_RE_k93_great_bedpe

less ../results/bedtools/revised_8c_enhancer/8cell_reprogrammed_enhancer_H3K27me3_for_GREAT_region_gene_table.txt \
   | sed -e 's, (,\t,g' | sed 's,)\,,\t,g' | sed 's,),,g' | sed 's,\t ,\t,g' \
   | awk '{print $1"\t"$2"\n"$1"\t"$4"\n"$1"\t"$6"\n"$1"\t"$8}' | awk '{if($2!="") print $0}' \
   | grep -v "^#" | sort -u | sort -k2,2 \
   | join -1 2 -2 5 - ../metadata/Homo_sapiens.GRCh38.97_pro3k_ucsc.bed \
   | sed 's, ,\t,g' | cut -f 2-5 \
   | sed -e 's,:,\t,g' -e 's,-,\t,g' > ../results/bedtools/revised_8c_enhancer/8cell_reprogrammed_enhancer_H3K27me3_for_GREAT_region_gene_table.bed
c8_RE_k273_great_bedpe=../results/bedtools/revised_8c_enhancer/8cell_reprogrammed_enhancer_H3K27me3_for_GREAT_region_gene_table.bedpe
cut -f 1-6 ../results/bedtools/revised_8c_enhancer/8cell_reprogrammed_enhancer_H3K27me3_for_GREAT_region_gene_table.bed | sort -k1,1 -k2,2n > $c8_RE_k273_great_bedpe

less ../results/bedtools/revised_8c_enhancer/8cell_primed_enhancer_for_GREAT_region_gene_table.txt \
   | sed -e 's, (,\t,g' | sed 's,)\,,\t,g' | sed 's,),,g' | sed 's,\t ,\t,g' \
   | awk '{print $1"\t"$2"\n"$1"\t"$4"\n"$1"\t"$6"\n"$1"\t"$8}' | awk '{if($2!="") print $0}' \
   | grep -v "^#" | sort -u | sort -k2,2 \
   | join -1 2 -2 5 - ../metadata/Homo_sapiens.GRCh38.97_pro3k_ucsc.bed \
   | sed 's, ,\t,g' | cut -f 2-5 \
   | sed -e 's,:,\t,g' -e 's,-,\t,g' > ../results/bedtools/revised_8c_enhancer/8cell_primed_enhancer_for_GREAT_region_gene_table.bed
c8_PE_great_bedpe=../results/bedtools/revised_8c_enhancer/8cell_primed_enhancer_for_GREAT_region_gene_table.bedpe
cut -f 1-6 ../results/bedtools/revised_8c_enhancer/8cell_primed_enhancer_for_GREAT_region_gene_table.bed | sort -k1,1 -k2,2n > $c8_PE_great_bedpe
# - process regions for aggregate plots
c8_RE_k93_bedpe=${od}/8cell_reprogrammed_enhancer_H3K9me3_interacted_with_genes.bedpe
c8_RE_k273_bedpe=${od}/8cell_reprogrammed_enhancer_H3K27me3_interacted_with_genes.bedpe
c8_PE_bepe=${od}/8cell_primed_enhancer_interacted_with_genes.bedpe
awk 'BEGIN{FS=OFS="\t"}{print $10,$11,$12,$13,$14,$15,".",$8}' ../results/bedtools/revised_8c_enhancer/genes_interacted_with_c8_RE_k93_enhancers_in_8cell_full.txt > ${c8_RE_k93_bedpe}
awk 'BEGIN{FS=OFS="\t"}{print $10,$11,$12,$13,$14,$15,".",$8}' ../results/bedtools/revised_8c_enhancer/genes_interacted_with_c8_RE_k273_enhancers_in_8cell_full.txt > ${c8_RE_k273_bedpe}
awk 'BEGIN{FS=OFS="\t"}{print $10,$11,$12,$13,$14,$15,".",$8}' ../results/bedtools/revised_8c_enhancer/genes_interacted_with_c8_PE_enhancers_in_8cell_full.txt > ${c8_PE_bepe}
sva_c8_RE_k93_bedpe=${od}/8cell_reprogrammed_SVA_enhancer_H3K9me3_interacted_with_genes.bedpe
sva_c8_RE_k273_bedpe=${od}/8cell_reprogrammed_SVA_enhancer_H3K27me3_interacted_with_genes.bedpe
sva_c8_PE_bepe=${od}/8cell_primed_SVA_enhancer_interacted_with_genes.bedpe
repeat_anno=/home/yhw/genome/ucsc/homo_sapiens/hg38/repeat/GRCh38_RepeatMasker_repeat_region.bed
less ${repeat_anno} | grep "SVA" > ./hs_SVA.bed
awk 'BEGIN{FS=OFS="\t"}{print $10,$11,$12,$1,$2,$3,$4,$5,$6,$7,$8,$9}' ../results/bedtools/revised_8c_enhancer/genes_interacted_with_c8_RE_k93_enhancers_in_8cell_full.txt \
   | bedtools intersect -a stdin -b ./hs_SVA.bed -wa | cut -f 4-12 > ${sva_c8_RE_k93_bedpe}
awk 'BEGIN{FS=OFS="\t"}{print $10,$11,$12,$1,$2,$3,$4,$5,$6,$7,$8,$9}' ../results/bedtools/revised_8c_enhancer/genes_interacted_with_c8_RE_k273_enhancers_in_8cell_full.txt \
   | bedtools intersect -a stdin -b ./hs_SVA.bed -wa | cut -f 4-12 > ${sva_c8_RE_k273_bedpe}
awk 'BEGIN{FS=OFS="\t"}{print $10,$11,$12,$1,$2,$3,$4,$5,$6,$7,$8,$9}' ../results/bedtools/revised_8c_enhancer/genes_interacted_with_c8_PE_enhancers_in_8cell_full.txt \
   | bedtools intersect -a stdin -b ./hs_SVA.bed -wa | cut -f 4-12 > ${sva_c8_PE_bepe}
rm ./hs_SVA.bed



### >>> Plotting signal on mouse ZGA genes (CMQ)
# - major ZGA gene
mm_zga_gene=/home/cmq/bioinfo/project-cmq/gse66582/results/R/Tables/zga/mm_zga_gene_genebody.bed
cat ${mm_zga_gene} | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,".",$5}' > metadata/mm_zga_gene_for_plot.bed
# - H3K9me3 bw
ocy_k93_fc_1=/home/cmq/bioinfo/project-cmq/gse98149_mouse_93/results/deeptools/bigwig/dedup_merged_bam_fc/SRR5479594_vs_SRR5479596.bw
ocy_k93_fc_2=/home/cmq/bioinfo/project-cmq/gse98149_mouse_93/results/deeptools/bigwig/dedup_merged_bam_fc/SRR5479595_vs_SRR5479596.bw
sp_k93_fc_1=/home/cmq/bioinfo/project-cmq/gse98149_mouse_93/results/deeptools/bigwig/dedup_merged_bam_fc/SRR5479597_vs_SRR5479596.bw
sp_k93_fc_2=/home/cmq/bioinfo/project-cmq/gse98149_mouse_93/results/deeptools/bigwig/dedup_merged_bam_fc/SRR5479598_vs_SRR5479596.bw
zyg_k93_fc=/home/cmq/bioinfo/project-cmq/gse98149_mouse_93/results/deeptools/bigwig/dedup_merged_bam_fc/SRR5479605_vs_SRR5479599.bw
c2_k93_fc=/home/cmq/bioinfo/project-cmq/gse98149_mouse_93/results/deeptools/bigwig/dedup_merged_bam_fc/SRR5479607_vs_SRR5479599.bw
# - H3K27me3 bw
c2_k273_fc=/home/cmq/bioinfo/project-cmq/gse73952_mouse_27/results/deeptools/bigwig/merged_bam_fc/SRR3208756_vs_SRR3208752.bw
ocy_k273_fc=/home/cmq/bioinfo/project-cmq/gse73952_mouse_27/results/deeptools/bigwig/merged_bam_fc/SRR3208749_vs_SRR3208744.bw
ps_k273_fc=/home/data/publicdata/GSE68507/analysis/mm_chipseq/results/deeptools/coverage/dedup_bam_fc/SRR1977633_trimmed_psort_dedup_vs_SRR1977634_trimmed_psort_dedup.bw
rs_k273_fc=/home/data/publicdata/GSE68507/analysis/mm_chipseq/results/deeptools/coverage/dedup_bam_fc/SRR1977637_trimmed_psort_dedup_vs_SRR1977638_trimmed_psort_dedup.bw
# - H3K27ac bw
c2_27ac_fc=/home/cmq/bioinfo/project-cmq/mm_27ac/results/deeptools/bigwig/merged_bam_fc/2cell_merge_vs_2cell_input_merge.bw
# - ATAC bw
sp_atac_1=/home/data/publicdata/GSE116857/analysis/atac/mouse/results/deeptools/coverage/SRR7504533_1_val_1_SRR7504533_2_val_2_psort_dedup.bw
sp_atac_2=/home/data/publicdata/GSE116857/analysis/atac/mouse/results/deeptools/coverage/SRR7504534_1_val_1_SRR7504534_2_val_2_psort_dedup.bw
gv_ocy_atac_1=/home/data/publicdata/GSE116857/analysis/atac/mouse/results/deeptools/coverage/SRR7504539_1_val_1_SRR7504539_2_val_2_psort_dedup.bw
gv_ocy_atac_2=/home/data/publicdata/GSE116857/analysis/atac/mouse/results/deeptools/coverage/SRR7504540_1_val_1_SRR7504540_2_val_2_psort_dedup.bw
m2_ocy_atac_1=/home/data/publicdata/GSE116857/analysis/atac/mouse/results/deeptools/coverage/SRR7504541_1_val_1_SRR7504541_2_val_2_psort_dedup.bw
m2_ocy_atac_2=/home/data/publicdata/GSE116857/analysis/atac/mouse/results/deeptools/coverage/SRR7504542_1_val_1_SRR7504542_2_val_2_psort_dedup.bw
zyg_atac_1=/home/data/publicdata/GSE136403/analysis/results/deeptools/coverage/dedup_bam/SRR10590670_1_val_1_SRR10590670_2_val_2_psort_dedup.bw
zyg_atac_2=/home/data/publicdata/GSE136403/analysis/results/deeptools/coverage/dedup_bam/SRR10590671_1_val_1_SRR10590671_2_val_2_psort_dedup.bw
c2_atac_1=/home/data/publicdata/GSE66581/analysis/results/deeptools/coverage/nfr_merge/2cell_rep1_merge.bw
c2_atac_2=/home/data/publicdata/GSE66581/analysis/results/deeptools/coverage/nfr_merge/2cell_rep2_merge.bw
e2_atac_1=/home/data/publicdata/GSE66581/analysis/results/deeptools/coverage/nfr_merge/e2cell_rep1_merge.bw
e2_atac_2=/home/data/publicdata/GSE66581/analysis/results/deeptools/coverage/nfr_merge/e2cell_rep2_merge.bw
# - H3K4me3 bw
gv_k43=/home/data/publicdata/GSE71434/analysis/chip/results/deeptools/coverage/SRR4014432_1_val_1_SRR4014432_2_val_2_psort_dedup.bw
m2o_k43_1=/home/data/publicdata/GSE71434/analysis/chip/results/deeptools/coverage/merge/m2o_rep1_merge.bw
m2o_k43_2=/home/data/publicdata/GSE71434/analysis/chip/results/deeptools/coverage/merge/m2o_rep2_merge.bw
sperm_k43_1=/home/data/publicdata/GSE71434/analysis/chip/results/deeptools/coverage/SRR3309386_1_val_1_SRR3309386_2_val_2_psort_dedup.bw
zyg_k43_1=/home/data/publicdata/GSE71434/analysis/chip/results/deeptools/coverage/merge/zygote_rep1_merge.bw
zyg_k43_2=/home/data/publicdata/GSE71434/analysis/chip/results/deeptools/coverage/merge/zygote_rep2_merge.bw
e2c_k43_1=/home/data/publicdata/GSE71434/analysis/chip/results/deeptools/coverage/SRR4014442_1_val_1_SRR4014442_2_val_2_psort_dedup.bw
e2c_k43_2=/home/data/publicdata/GSE71434/analysis/chip/results/deeptools/coverage/SRR4014443_1_val_1_SRR4014443_2_val_2_psort_dedup.bw
l2c_k43_1=/home/data/publicdata/GSE71434/analysis/chip/results/deeptools/coverage/SRR4014444_1_val_1_SRR4014444_2_val_2_psort_dedup.bw
l2c_k43_2=/home/data/publicdata/GSE71434/analysis/chip/results/deeptools/coverage/SRR4014445_1_val_1_SRR4014445_2_val_2_psort_dedup.bw
# RNA PolII bw
m2_ooc_pol2_1=/home/data/publicdata/GSE135457/analysis/gse135457_2/results/deeptools/coverage/for_k9/MII_rep1_cp.bw
m2_ooc_pol2_2=/home/data/publicdata/GSE135457/analysis/gse135457_2/results/deeptools/coverage/for_k9/MII_rep2_cp.bw
pn3_pol2_1=/home/data/publicdata/GSE135457/analysis/gse135457_2/results/deeptools/coverage/for_k9/PN3_rep1_cp.bw
pn3_pol2_2=/home/data/publicdata/GSE135457/analysis/gse135457_2/results/deeptools/coverage/for_k9/PN3_rep2_cp.bw
pn5_pol2_1=/home/data/publicdata/GSE135457/analysis/gse135457_2/results/deeptools/coverage/for_k9/PN5_Pol2_rep1_merge.bw
pn5_pol2_2=/home/data/publicdata/GSE135457/analysis/gse135457_2/results/deeptools/coverage/for_k9/PN5_Pol2_rep2_merge.bw
e2c_pol2_1=/home/data/publicdata/GSE135457/analysis/gse135457_2/results/deeptools/coverage/for_k9/2cell_early_Pol2_rep1_merge.bw
e2c_pol2_2=/home/data/publicdata/GSE135457/analysis/gse135457_2/results/deeptools/coverage/for_k9/2cell_early_Pol2_rep2_merge.bw
l2c_pol2_1=/home/data/publicdata/GSE135457/analysis/gse135457_2/results/deeptools/coverage/for_k9/2cell_late_Pol2_rep1_merge.bw
l2c_pol2_2=/home/data/publicdata/GSE135457/analysis/gse135457_2/results/deeptools/coverage/for_k9/2cell_late_Pol2_rep2_merge.bw
l2c_pol2_3=/home/data/publicdata/GSE135457/analysis/gse135457_2/results/deeptools/coverage/for_k9/2cell_late_Pol2_rep3_merge.bw
# - create output dir
mkdir -p results/deeptools/matrix/revised_mm_zga_gene
mkdir -p results/deeptools/heatmap/revised_mm_zga_gene
# - plotting
computeMatrix scale-regions -S ${ocy_k93_fc_1} ${ocy_k93_fc_2} ${sp_k93_fc_1} ${sp_k93_fc_2} ${zyg_k93_fc} ${c2_k93_fc} \
                               ${ocy_k273_fc} ${ps_k273_fc} ${rs_k273_fc} ${c2_k273_fc} ${c2_27ac_fc} \
                               ${gv_k43} ${m2o_k43_1} ${m2o_k43_2} ${sperm_k43_1} ${zyg_k43_1} ${zyg_k43_2} ${e2c_k43_1} ${e2c_k43_2} ${l2c_k43_1} ${l2c_k43_2}\
                               ${gv_ocy_atac_1} ${gv_ocy_atac_2} ${m2_ocy_atac_1} ${m2_ocy_atac_2} ${sp_atac_1} ${sp_atac_2} \
                               ${zyg_atac_1} ${zyg_atac_2} ${e2_atac_1} ${e2_atac_2} ${c2_atac_1} ${c2_atac_2} \
                               ${m2_ooc_pol2_1} ${m2_ooc_pol2_2} ${pn3_pol2_1} ${pn3_pol2_2} ${pn5_pol2_1} ${pn5_pol2_2} ${e2c_pol2_1} ${e2c_pol2_2} ${l2c_pol2_1} ${l2c_pol2_2} ${l2c_pol2_3} \
                            -R metadata/mm_zga_gene_for_plot.bed --averageTypeBins mean --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
                            --numberOfProcessors 30 --missingDataAsZero --metagene \
                            -o results/deeptools/matrix/revised_mm_zga_gene/gzmatrix.gz
plotHeatmap -m results/deeptools/matrix/revised_mm_zga_gene/gzmatrix.gz \
            -o results/deeptools/heatmap/revised_mm_zga_gene/mm_zga_genes_H3K9me3_H3K27me3_H3K27ac_H3K4me3_ATAC_RNAPol2_heatmap_-2-4.pdf \
            --averageTypeSummaryPlot mean  --heatmapHeight 6 --heatmapWidth 2 --whatToShow "plot, heatmap and colorbar" \
            --plotType lines --colorList cyan,white,red --zMin -2 --zMax 4 --yMin -2.5 --yMax 2.5 --legendLocation upper-left --interpolationMethod auto --missingDataColor 1 \
            --samplesLabel "ocy_k93_1" "ocy_k93_2" "sperm_k93_1" "sperm_k93_2" "zyg_k93" "2c_k93" \
                           "ocy_k273" "ps_k273" "rs_k273" "2c_k273" "2c_k27ac" \
                           "gv_k43" "m2o_k43_1" "m2o_k43_2" "sperm_k43_1" "zyg_k43_1" "zyg_k43_2" "e2c_k43_1" "e2c_k43_2" "l2c_k43_1" "l2c_k43_2" \
                           "gv_atac_1" "gv_atac_2" "m2o_atac_1" "m2o_atac_2" "sperm_atac_1" "sperm_atac_2" \
                           "zyg_atac_1" "zyg_atac_2" "e2c_atac_1" "e2c_atac_2" "2c_atac_1" "2c_atac_2" \
                           "m2_ooc_pol2_1" "m2_ooc_pol2_2" "pn3_pol2_1" "pn3_pol2_2" "pn5_pol2_1" "pn5_pol2_2" \
                           "e2c_pol2_1" "e2c_pol2_2" "l2c_pol2_1" "l2c_pol2_2" "l2c_pol2_3" \
            --regionsLabel "Mouse Major ZGA Genes"
for max in `seq 5 5 100`
do
    plotHeatmap -m results/deeptools/matrix/revised_mm_zga_gene/gzmatrix.gz \
                -o results/deeptools/heatmap/revised_mm_zga_gene/mm_zga_genes_H3K9me3_H3K27me3_H3K27ac_H3K4me3_ATAC_RNAPol2_heatmap_0-${max}.pdf \
                --averageTypeSummaryPlot mean  --heatmapHeight 6 --heatmapWidth 2 --whatToShow "plot, heatmap and colorbar" \
                --plotType lines --colorMap Blues --zMin 0 --zMax ${max} --legendLocation upper-left --interpolationMethod auto --missingDataColor 1 \
                --samplesLabel "ocy_k93_1" "ocy_k93_2" "sperm_k93_1" "sperm_k93_2" "zyg_k93" "2c_k93" \
                               "ocy_k273" "ps_k273" "rs_k273" "2c_k273" "2c_k27ac" \
                               "gv_k43" "m2o_k43_1" "m2o_k43_2" "sperm_k43_1" "zyg_k43_1" "zyg_k43_2" "e2c_k43_1" "e2c_k43_2" "l2c_k43_1" "l2c_k43_2" \
                               "gv_atac_1" "gv_atac_2" "m2o_atac_1" "m2o_atac_2" "sperm_atac_1" "sperm_atac_2" \
                               "zyg_atac_1" "zyg_atac_2" "e2c_atac_1" "e2c_atac_2" "2c_atac_1" "2c_atac_2" \
                               "m2_ooc_pol2_1" "m2_ooc_pol2_2" "pn3_pol2_1" "pn3_pol2_2" "pn5_pol2_1" "pn5_pol2_2" \
                               "e2c_pol2_1" "e2c_pol2_2" "l2c_pol2_1" "l2c_pol2_2" "l2c_pol2_3" \
                --regionsLabel "Mouse Major ZGA Genes"
done



### >>> Plotting signal on human ZGA genes (CMQ)
# - human zga gene
hs_zga_gene=/home/cmq/bioinfo/project-cmq/rna_seq/results/R/Tables/zga/hs_zga_gene_genebody.bed
cat ${hs_zga_gene} | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,".",$5}' > metadata/hs_zga_gene_for_plot.bed
# - H3K9me bigwig files
sperm_k93=/home/cmq/bioinfo/project-cmq/embryo_93/results/deeptools/bigwig/merged_bam/H3K9me3_rep1_trimmed_psort_dedup.bw
sperm_k93=/home/cmq/bioinfo/project-cmq/embryo_93/results/deeptools/bigwig/merged_bam/H3K9me3_rep2_trimmed_psort_dedup.bw
c4_k93_fc=/home/cmq/bioinfo/project-cmq/embryo_93/results/deeptools/bigwig/merged_bam_ratio/h4C93-Ip-1_vs_h4C93-Input-1.bw
c8_k93_fc=/home/cmq/bioinfo/project-cmq/embryo_93/results/deeptools/bigwig/merged_bam_ratio/h8C93-Ip-3_vs_h8C93-Input-2.bw
# - H3K27me3 bigwig files
oocyte_k273=/home/cmq/bioinfo/project-cmq/gse124718/results/deeptools/bigwig/merge_replicate_multi/SRR8402604_1_val_1_SRR8402604_2_val_2_psorted_dedup_SRR8402605_1_val_1_SRR8402605_2_val_2_psorted_dedup_merged.bw
ps_k273_fc_1=/home/data/publicdata/GSE68507/analysis/human/results/deeptools/coverage/foldchange/SRR1977543_human_1_pachytene_spermatocytes_H3K27me3_ChIP-seq_trimmed_psort_dedup_vs_SRR1977545_human_1_pachytene_spermatocyte_input_for_ChIP_trimmed_psort_dedup.bw
ps_k273_fc_2=/home/data/publicdata/GSE68507/analysis/human/results/deeptools/coverage/foldchange/SRR1977600_human_2_pachytene_spermatocyte_H3K27me3_ChIP-seq_trimmed_psort_dedup_vs_SRR1977602_human_2_pachytene_spermatocyte_input_for_ChIP_trimmed_psort_dedup.bw
ps_k273_fc_3=/home/data/publicdata/GSE68507/analysis/human/results/deeptools/coverage/foldchange/SRR1977613_human_3_pachytene_spermatocyte_H3K27me3_ChIP-seq_trimmed_psort_dedup_vs_SRR1977614_human_3_pachytene_spermatocyte_input_for_ChIP_trimmed_psort_dedup.bw
rs_k273_fc_1=/home/data/publicdata/GSE68507/analysis/human/results/deeptools/coverage/foldchange/SRR1977584_human_1_round_spermatid_H3K27me3_ChIP-seq_trimmed_psort_dedup_vs_SRR1977586_human_1_round_spermatid_input_for_ChIP_trimmed_psort_dedup.bw
rs_k273_fc_2=/home/data/publicdata/GSE68507/analysis/human/results/deeptools/coverage/foldchange/SRR1977606_human_2_round_spermatid_H3K27me3_ChIP-seq_trimmed_psort_dedup_vs_SRR1977608_human_2_round_spermatid_input_for_ChIP_trimmed_psort_dedup.bw
rs_k273_fc_3=/home/data/publicdata/GSE68507/analysis/human/results/deeptools/coverage/foldchange/SRR1977617_human_3_round_spermatid_H3K27me3_ChIP-seq_trimmed_psort_dedup_vs_SRR1977614_human_3_pachytene_spermatocyte_input_for_ChIP_trimmed_psort_dedup.bw
c2_k273=/home/cmq/bioinfo/project-cmq/gse124718/results/deeptools/bigwig/merge_replicate_multi/2cell_multi_replicate_merged.bw
c4_k273=/home/cmq/bioinfo/project-cmq/gse124718/results/deeptools/bigwig/merge_replicate_multi/4cell_multi_replicate_merged.bw
c8_k273=/home/cmq/bioinfo/project-cmq/gse124718/results/deeptools/bigwig/merge_replicate_multi/8cell_multi_replicate_merged.bw
# - H3K27ac bigwig files
c8_k27ac_1=/home/yhw/bioinfo/project-temp/wjc/gse124718/results/deeptools/coverage/SRR9131735_1_val_1_SRR9131735_2_val_2_psort_dedup.bw
c8_k27ac_2=/home/yhw/bioinfo/project-temp/wjc/gse124718/results/deeptools/coverage/SRR9131736_1_val_1_SRR9131736_2_val_2_psort_dedup.bw
# - H3K4me3 bigwig files
gv_oocyte_k43=/home/cmq/bioinfo/project-cmq/gse124718_human_43/results/deeptools/bigwig/merged_bam/SRR8402601_1_val_1_SRR8402601_2_val_2_psorted_dedup_SRR8402602_1_val_1_SRR8402602_2_val_2_psorted_dedup_merged.bw
m2_oocyte_k43=/home/cmq/bioinfo/project-cmq/gse124718_human_43/results/deeptools/bigwig/raw_bam/SRR8402606_1_val_1_SRR8402606_2_val_2_psorted_dedup.bw
ps_k43_fc_1=/home/data/publicdata/GSE68507/analysis/human/results/deeptools/coverage/foldchange/SRR1977542_human_1_pachytene_spermatocyte_H3K4me3_ChIP-seq_trimmed_psort_dedup_vs_SRR1977545_human_1_pachytene_spermatocyte_input_for_ChIP_trimmed_psort_dedup.bw
ps_k43_fc_2=/home/data/publicdata/GSE68507/analysis/human/results/deeptools/coverage/foldchange/SRR1977598_human_2_pachytene_spermatocyte_H3K4me3_ChIP-seq_trimmed_psort_dedup_vs_SRR1977602_human_2_pachytene_spermatocyte_input_for_ChIP_trimmed_psort_dedup.bw
ps_k43_fc_3=/home/data/publicdata/GSE68507/analysis/human/results/deeptools/coverage/foldchange/SRR1977612_human_3_pachytene_spermatocyte_H3K4me3_ChIP-seq_trimmed_psort_dedup_vs_SRR1977614_human_3_pachytene_spermatocyte_input_for_ChIP_trimmed_psort_dedup.bw
rs_k43_fc_1=/home/data/publicdata/GSE68507/analysis/human/results/deeptools/coverage/foldchange/SRR1977583_human_1_round_spermatid_H3K4me3_ChIP-seq_trimmed_psort_dedup_vs_SRR1977586_human_1_round_spermatid_input_for_ChIP_trimmed_psort_dedup.bw
rs_k43_fc_2=/home/data/publicdata/GSE68507/analysis/human/results/deeptools/coverage/foldchange/SRR1977605_human_2_round_spermatid_H3K4me3_ChIP-seq_trimmed_psort_dedup_vs_SRR1977608_human_2_round_spermatid_input_for_ChIP_trimmed_psort_dedup.bw
rs_k43_fc_3=/home/data/publicdata/GSE68507/analysis/human/results/deeptools/coverage/foldchange/SRR1977616_human_3_round_spermatid_H3K4me3_ChIP-seq_trimmed_psort_dedup_vs_SRR1977614_human_3_pachytene_spermatocyte_input_for_ChIP_trimmed_psort_dedup.bw
c4_k43=/home/cmq/bioinfo/project-cmq/gse124718_human_43/results/deeptools/bigwig/merged_bam/SRR8402611_1_val_1_SRR8402611_2_val_2_psorted_dedup_SRR8402612_trimmed_psorted_dedup_merged.bw
c8_k43=/home/cmq/bioinfo/project-cmq/gse124718_human_43/results/deeptools/bigwig/merged_bam/SRR8402624_1_val_SRR8402625_1_val_SRR8402627_1_val_SRR8402628_1_valmerged.bw
# - ATAC bigwig files
gv_oocyte_atac_1=/home/data/publicdata/GSE124718/analysis/atac/results/deeptools/coverage/SRR9131727_1_val_1_SRR9131727_2_val_2_psort_dedup.bw
gv_oocyte_atac_2=/home/data/publicdata/GSE124718/analysis/atac/results/deeptools/coverage/SRR9131728_1_val_1_SRR9131728_2_val_2_psort_dedup.bw
sperm_atac_1=/home/data/publicdata/GSE116857/analysis/atac/human/results/deeptools/coverage/SRR7504537_1_val_1_SRR7504537_2_val_2_psort_dedup.bw
sperm_atac_2=/home/data/publicdata/GSE116857/analysis/atac/human/results/deeptools/coverage/SRR7504538_1_val_1_SRR7504538_2_val_2_psort_dedup.bw
zy_atac=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/deeptools/bigwig/dedup_merge_bam/rename/zygote_3pn_dedup.bw
c2_atac_1=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/deeptools/bigwig/dedup_merge_bam/rename/c2_3pn_1_dedup_merged.bw
c2_atac_2=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/deeptools/bigwig/dedup_merge_bam/rename/c2_3pn_2_dedup_merged.bw
c4_atac_1=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/deeptools/bigwig/dedup_merge_bam/SRR5837324_1_val_1_SRR5837324_2_val_2_psorted_dedup_SRR5837325_1_val_1_SRR5837325_2_val_2_psorted_dedup_merged.bw
c4_atac_2=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/deeptools/bigwig/dedup_merge_bam/SRR5837326_1_val_1_SRR5837326_2_val_2_psorted_dedup_SRR5837327_1_val_1_SRR5837327_2_val_2_psorted_dedup_SRR5837328_1_val_1_SRR5837328_2_val_2_psorted_dedup_merged.bw
c8_atac_1=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/deeptools/bigwig/dedup_merge_bam/SRR5837333_1_val_1_SRR5837333_2_val_2_psorted_dedup_SRR5837334_1_val_1_SRR5837334_2_val_2_psorted_dedup_merged.bw
c8_atac_2=/home/cmq/bioinfo/project-ymz/atac/gse101571/results/deeptools/bigwig/dedup_merge_bam/SRR5837335_1_val_1_SRR5837335_2_val_2_psorted_dedup_SRR5837336_1_val_1_SRR5837336_2_val_2_psorted_dedup_merged.bw
# - plotting
mkdir -p results/deeptools/matrix/revised_hs_zga_gene
mkdir -p results/deeptools/heatmap/revised_hs_zga_gene
#<<'DONE'
computeMatrix scale-regions -S ${sperm_k93} ${sperm_k93} ${c4_k93_fc} ${c8_k93_fc} \
                               ${oocyte_k273} ${ps_k273_fc_1} ${ps_k273_fc_2} ${ps_k273_fc_3} ${rs_k273_fc_1} ${rs_k273_fc_2} ${rs_k273_fc_3} ${c2_k273} ${c4_k273} ${c8_k273} \
                               ${c8_k27ac_1} ${c8_k27ac_2} ${gv_oocyte_k43} ${m2_oocyte_k43} ${ps_k43_fc_1} ${ps_k43_fc_2} ${ps_k43_fc_3} ${rs_k43_fc_1} ${rs_k43_fc_2} ${rs_k43_fc_3} ${c4_k43} ${c8_k43} \
                               ${gv_oocyte_atac_1} ${gv_oocyte_atac_2} ${sperm_atac_1} ${sperm_atac_2} ${zy_atac} ${c2_atac_1} ${c2_atac_2} ${c4_atac_1} ${c4_atac_2} ${c8_atac_1} ${c8_atac_2} \
                            -R metadata/hs_zga_gene_for_plot.bed --averageTypeBins mean --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
                            --numberOfProcessors 20 --missingDataAsZero --binSize 10 --metagene \
                            -o results/deeptools/matrix/revised_hs_zga_gene/gzmatrix.gz
# ATAC
for max in `seq 20 5 45`
do
    plotHeatmap -m results/deeptools/matrix/revised_hs_zga_gene/gzmatrix.gz \
                -o results/deeptools/heatmap/revised_hs_zga_gene/hs_zga_gene_H3K9me3_H3K27me3_H3K27ac_H3K4me3_and_ATAC_signal_0-${max}_ATAC_heatmap.pdf \
                --averageTypeSummaryPlot mean  --heatmapHeight 6 --heatmapWidth 2 --whatToShow "plot, heatmap and colorbar" \
                --plotType lines --colorMap Blues --zMin 0 --zMax ${max} \
                --legendLocation upper-left --interpolationMethod auto --missingDataColor 1 \
                --samplesLabel "sperm_k93" "sperm_k93" "c4_k93_fc" "c8_k93_fc" \
                               "ooc_k273" "ps_k273_fc_1" "ps_k273_fc_2" "ps_k273_fc_3" "rs_k273_fc_1" "rs_k273_fc_2" "rs_k273_fc_3" "c2_k273" "c4_k273" "c8_k273" \
                               "c8_k27ac_1" "c8_k27ac_2" "gv_ooc_k43" "m2_ooc_k43" "ps_k43_fc_1" "ps_k43_fc_2" "ps_k43_fc_3" "rs_k43_fc_1" "rs_k43_fc_2" "rs_k43_fc_3" "c4_k43" "c8_k43" \
                               "gv_ooc_atac_1" "gv_ooc_atac_2" "sperm_atac_1" "sperm_atac2" "zy_atac" "c2_atac_1" "c2_atac_2" "c4_atac_1" "c4_atac_2" "c8_atac_1" "c8_atac_2" \
                --regionsLabel "human major zga gene"
done
# H3K27me3
for max in `seq 5 5 15`
do
    plotHeatmap -m results/deeptools/matrix/revised_hs_zga_gene/gzmatrix.gz \
                -o results/deeptools/heatmap/revised_hs_zga_gene/hs_zga_gene_H3K9me3_H3K27me3_H3K27ac_H3K4me3_and_ATAC_signal_0-${max}_H3K27me3_and_H3K9me3_sperm_heatmap.pdf \
                --averageTypeSummaryPlot mean  --heatmapHeight 6 --heatmapWidth 2 --whatToShow "plot, heatmap and colorbar" \
                --plotType lines --colorMap Blues --zMin 0 --zMax ${max} --yMin 0 --yMax 10 \
                --legendLocation upper-left --interpolationMethod auto --missingDataColor 1 \
                --samplesLabel "sperm_k93" "sperm_k93" "c4_k93_fc" "c8_k93_fc" \
                               "ooc_k273" "ps_k273_fc_1" "ps_k273_fc_2" "ps_k273_fc_3" "rs_k273_fc_1" "rs_k273_fc_2" "rs_k273_fc_3" "c2_k273" "c4_k273" "c8_k273" \
                               "c8_k27ac_1" "c8_k27ac_2" "gv_ooc_k43" "m2_ooc_k43" "ps_k43_fc_1" "ps_k43_fc_2" "ps_k43_fc_3" "rs_k43_fc_1" "rs_k43_fc_2" "rs_k43_fc_3" "c4_k43" "c8_k43" \
                               "gv_ooc_atac_1" "gv_ooc_atac_2" "sperm_atac_1" "sperm_atac2" "zy_atac" "c2_atac_1" "c2_atac_2" "c4_atac_1" "c4_atac_2" "c8_atac_1" "c8_atac_2" \
                --regionsLabel "human major zga gene"
done
# H3K27ac 
for max in `seq 40 5 60`
do
    plotHeatmap -m results/deeptools/matrix/revised_hs_zga_gene/gzmatrix.gz \
                -o results/deeptools/heatmap/revised_hs_zga_gene/hs_zga_gene_H3K9me3_H3K27me3_H3K27ac_H3K4me3_and_ATAC_signal_0-${max}_H3K27ac_heatmap.pdf \
                --averageTypeSummaryPlot mean  --heatmapHeight 6 --heatmapWidth 2 --whatToShow "plot, heatmap and colorbar" \
                --plotType lines --colorMap Greens --zMin 0 --zMax ${max} --yMin 0 --yMax 50 \
                --legendLocation upper-left --interpolationMethod auto --missingDataColor 1 \
                --samplesLabel "sperm_k93" "sperm_k93" "c4_k93_fc" "c8_k93_fc" \
                               "ooc_k273" "ps_k273_fc_1" "ps_k273_fc_2" "ps_k273_fc_3" "rs_k273_fc_1" "rs_k273_fc_2" "rs_k273_fc_3" "c2_k273" "c4_k273" "c8_k273" \
                               "c8_k27ac_1" "c8_k27ac_2" "gv_ooc_k43" "m2_ooc_k43" "ps_k43_fc_1" "ps_k43_fc_2" "ps_k43_fc_3" "rs_k43_fc_1" "rs_k43_fc_2" "rs_k43_fc_3" "c4_k43" "c8_k43" \
                               "gv_ooc_atac_1" "gv_ooc_atac_2" "sperm_atac_1" "sperm_atac2" "zy_atac" "c2_atac_1" "c2_atac_2" "c4_atac_1" "c4_atac_2" "c8_atac_1" "c8_atac_2" \
                --regionsLabel "human major zga gene"
done
# H3K4me3
for max in `seq 50 5 100`
do
    plotHeatmap -m results/deeptools/matrix/revised_hs_zga_gene/gzmatrix.gz \
                -o results/deeptools/heatmap/revised_hs_zga_gene/hs_zga_gene_H3K9me3_H3K27me3_H3K27ac_H3K4me3_and_ATAC_signal_0-${max}_H3K4me3_heatmap.pdf \
                --averageTypeSummaryPlot mean  --heatmapHeight 6 --heatmapWidth 2 --whatToShow "plot, heatmap and colorbar" \
                --plotType lines --colorMap Oranges --zMin 0 --zMax ${max} --yMin 0 --yMax 90\
                --legendLocation upper-left --interpolationMethod auto --missingDataColor 1 \
                --samplesLabel "sperm_k93" "sperm_k93" "c4_k93_fc" "c8_k93_fc" \
                               "ooc_k273" "ps_k273_fc_1" "ps_k273_fc_2" "ps_k273_fc_3" "rs_k273_fc_1" "rs_k273_fc_2" "rs_k273_fc_3" "c2_k273" "c4_k273" "c8_k273" \
                               "c8_k27ac_1" "c8_k27ac_2" "gv_ooc_k43" "m2_ooc_k43" "ps_k43_fc_1" "ps_k43_fc_2" "ps_k43_fc_3" "rs_k43_fc_1" "rs_k43_fc_2" "rs_k43_fc_3" "c4_k43" "c8_k43" \
                               "gv_ooc_atac_1" "gv_ooc_atac_2" "sperm_atac_1" "sperm_atac2" "zy_atac" "c2_atac_1" "c2_atac_2" "c4_atac_1" "c4_atac_2" "c8_atac_1" "c8_atac_2" \
                --regionsLabel "human major zga gene"
done
# foldchange bw
for max in `seq 1 1 5`
do
    plotHeatmap     -m results/deeptools/matrix/revised_hs_zga_gene/gzmatrix.gz \
                    -o results/deeptools/heatmap/revised_hs_zga_gene/hs_zga_gene_H3K9me3_H3K27me3_H3K27ac_H3K4me3_and_ATAC_signal_-1-${max}_foldchange_heatmap.pdf \
                    --averageTypeSummaryPlot mean  --heatmapHeight 6 --heatmapWidth 2 --whatToShow "plot, heatmap and colorbar" \
                    --plotType lines --colorList cyan,white,red --zMin -1 --zMax ${max} --yMin -1.5 --yMax 2.5 \
                    --legendLocation upper-left --interpolationMethod auto --missingDataColor 1 \
                    --samplesLabel "sperm_k93" "sperm_k93" "c4_k93_fc" "c8_k93_fc" \
                                   "ooc_k273" "ps_k273_fc_1" "ps_k273_fc_2" "ps_k273_fc_3" "rs_k273_fc_1" "rs_k273_fc_2" "rs_k273_fc_3" "c2_k273" "c4_k273" "c8_k273" \
                                   "c8_k27ac_1" "c8_k27ac_2" "gv_ooc_k43" "m2_ooc_k43" "ps_k43_fc_1" "ps_k43_fc_2" "ps_k43_fc_3" "rs_k43_fc_1" "rs_k43_fc_2" "rs_k43_fc_3" "c4_k43" "c8_k43" \
                                   "gv_ooc_atac_1" "gv_ooc_atac_2" "sperm_atac_1" "sperm_atac2" "zy_atac" "c2_atac_1" "c2_atac_2" "c4_atac_1" "c4_atac_2" "c8_atac_1" "c8_atac_2" \
                    --regionsLabel "human major zga gene"
done



### >>> Plot Hic aggregate (YHW)
c8_PE_bedpe=/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer/8cell_primed_enhancer_interacted_with_genes.bedpe
c8_RE_k9_bedpe=/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer/8cell_reprogrammed_enhancer_H3K9me3_interacted_with_genes.bedpe
c8_RE_k27_bedpe=/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer/8cell_reprogrammed_enhancer_H3K27me3_interacted_with_genes.bedpe
c8_RE_bedpe=./results/8cell_reprogrammed_enhancer_interacted_with_genes.bedpe
c8_PE_great_bedpe=/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer/8cell_primed_enhancer_for_GREAT_region_gene_table.bedpe
c8_RE_k9_great_bedpe=/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer/8cell_reprogrammed_enhancer_H3K9me3_for_GREAT_region_gene_table.bedpe
c8_RE_k27_great_bedpe=/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer/8cell_reprogrammed_enhancer_H3K27me3_for_GREAT_region_gene_table.bedpe
sva_c8_PE_bepe=/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer/8cell_primed_SVA_enhancer_interacted_with_genes.bedpe
sva_c8_RE_k273_bedpe=/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer/8cell_reprogrammed_SVA_enhancer_H3K27me3_interacted_with_genes.bedpe
sva_c8_RE_k93_bedpe=/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer/8cell_reprogrammed_SVA_enhancer_H3K9me3_interacted_with_genes.bedpe
sva_c8_bedpe=./results/8cell_sva_enhancer_interacted_with_genes.bedpe
cat ${sva_c8_RE_k273_bedpe} ${sva_c8_RE_k93_bedpe} > ${sva_c8_bedpe}
cat ${c8_RE_k9_bedpe} ${c8_RE_k27_bedpe} > ${c8_RE_bedpe}
for bedpe in sva_c8_bedpe c8_RE_bedpe c8_PE_bedpe c8_RE_k9_bedpe c8_RE_k27_bedpe 
do
    enhancer=./results/enhancer_${bedpe}.bed
    cut -f 1-3 ${!bedpe} > ${enhancer}
    gene=./results/gene_${bedpe}.bed
    cut -f 4-6 ${!bedpe} > ${gene}
    # 8cell rep1
    hicAggregateContacts --matrix ./results/hicexplorer/matrix/correct/8cell_rep1_f1_bam_merge_nsort_8cell_rep1_r2_bam_merge_nsort_40kbnor_ICE.h5 \
                         --BED ${enhancer} --BED2 ${gene} --outFileContactPairs test_tmp \
                         --outFileName 8c_rep1_aggregate_Contacts_${bedpe}_0p5to1p5.pdf --vMin 0.5 --vMax 1.5 \
                         --range 1000:20000000 --numberOfBins 100 \
                         --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
                         --avgType mean --transform obs/exp --howToCluster center --colorMap bwr 
    hicAggregateContacts --matrix ./results/hicexplorer/matrix/correct/8cell_rep1_f1_bam_merge_nsort_8cell_rep1_r2_bam_merge_nsort_40kbnor_ICE.h5 \
                         --BED ${enhancer} --BED2 ${gene} \
                         --outFileName 8c_rep1_aggregate_Contacts_${bedpe}_0p5to1p75.pdf --vMin 0.5 --vMax 1.75 \
                         --range 1000:20000000 --numberOfBins 100 \
                         --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
                         --avgType mean --transform obs/exp --howToCluster center --colorMap bwr &
    hicAggregateContacts --matrix ./results/hicexplorer/matrix/correct/8cell_rep1_f1_bam_merge_nsort_8cell_rep1_r2_bam_merge_nsort_40kbnor_ICE.h5 \
                         --BED ${enhancer} --BED2 ${gene} \
                         --outFileName 8c_rep1_aggregate_Contacts_${bedpe}_0p5to2p0.pdf --vMin 0.5 --vMax 2.0 \
                         --range 1000:20000000 --numberOfBins 100 \
                         --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
                         --avgType mean --transform obs/exp --howToCluster center --colorMap bwr
    hicAggregateContacts --matrix ./results/hicexplorer/matrix/correct/8cell_rep1_f1_bam_merge_nsort_8cell_rep1_r2_bam_merge_nsort_40kbnor_ICE.h5 \
                         --BED ${enhancer} --BED2 ${gene} \
                         --outFileName 8c_rep1_aggregate_Contacts_${bedpe}_0p5to2p25.pdf --vMin 0.5 --vMax 2.25 \
                         --range 1000:20000000 --numberOfBins 100 \
                         --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
                         --avgType mean --transform obs/exp --howToCluster center --colorMap bwr &
    hicAggregateContacts --matrix ./results/hicexplorer/matrix/correct/8cell_rep1_f1_bam_merge_nsort_8cell_rep1_r2_bam_merge_nsort_40kbnor_ICE.h5 \
                         --BED ${enhancer} --BED2 ${gene} \
                         --outFileName 8c_rep1_aggregate_Contacts_${bedpe}_0p5to2p5.pdf --vMin 0.5 --vMax 2.5 \
                         --range 1000:20000000 --numberOfBins 100 \
                         --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
                         --avgType mean --transform obs/exp --howToCluster center --colorMap bwr
    # 8cell rep2
    hicAggregateContacts --matrix ./results/hicexplorer/matrix/correct/8cell_rep2_f1_bam_merge_nsort_8cell_rep2_r2_bam_merge_nsort_40kbnor_ICE.h5 \
                         --BED ${enhancer} --BED2 ${gene} \
                         --outFileName 8c_rep2_aggregate_Contacts_${bedpe}_0p5to1p5.pdf --vMin 0.5 --vMax 1.5 \
                         --range 1000:20000000 --numberOfBins 100 \
                         --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
                         --avgType mean --transform obs/exp --howToCluster center --colorMap bwr &
    hicAggregateContacts --matrix ./results/hicexplorer/matrix/correct/8cell_rep2_f1_bam_merge_nsort_8cell_rep2_r2_bam_merge_nsort_40kbnor_ICE.h5 \
                         --BED ${enhancer} --BED2 ${gene} \
                         --outFileName 8c_rep2_aggregate_Contacts_${bedpe}_0p5to1p75.pdf --vMin 0.5 --vMax 1.75 \
                         --range 1000:20000000 --numberOfBins 100 \
                         --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
                         --avgType mean --transform obs/exp --howToCluster center --colorMap bwr &
    hicAggregateContacts --matrix ./results/hicexplorer/matrix/correct/8cell_rep2_f1_bam_merge_nsort_8cell_rep2_r2_bam_merge_nsort_40kbnor_ICE.h5 \
                         --BED ${enhancer} --BED2 ${gene} \
                         --outFileName 8c_rep2_aggregate_Contacts_${bedpe}_0p5to2p0.pdf --vMin 0.5 --vMax 2.0 \
                         -range 1000:20000000 --numberOfBins 100 \
                         --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
                         --avgType mean --transform obs/exp --howToCluster center --colorMap bwr
    hicAggregateContacts --matrix ./results/hicexplorer/matrix/correct/8cell_rep2_f1_bam_merge_nsort_8cell_rep2_r2_bam_merge_nsort_40kbnor_ICE.h5 \
                         --BED ${enhancer} --BED2 ${gene} \
                         --outFileName 8c_rep2_aggregate_Contacts_${bedpe}_0p5to2p25.pdf --vMin 0.5 --vMax 2.25 \
                         --range 1000:20000000 --numberOfBins 100 \
                         --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
                         --avgType mean --transform obs/exp --howToCluster center --colorMap bwr &
    hicAggregateContacts --matrix ./results/hicexplorer/matrix/correct/8cell_rep2_f1_bam_merge_nsort_8cell_rep2_r2_bam_merge_nsort_40kbnor_ICE.h5 \
                         --BED ${enhancer} --BED2 ${gene} \
                         --outFileName 8c_rep2_aggregate_Contacts_${bedpe}_0p5to2p5.pdf --vMin 0.5 --vMax 2.5 \
                         --range 1000:20000000 --numberOfBins 100 \
                         --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
                         --avgType mean --transform obs/exp --howToCluster center --colorMap bwr
    # 2 cell rep1
    hicAggregateContacts --matrix results/hicexplorer/matrix/correct_all_40kb/2cell_rep1_f1_2cell_rep1_r2_40kb_nor_ICE.h5 \
                         --BED ${enhancer} --BED2 ${gene} \
                         --outFileName 2c_rep1_aggregate_Contacts_${bedpe}_0p5to1p5.pdf --vMin 0.5 --vMax 1.5 \
                         --range 1000:20000000 --numberOfBins 100 \
                         --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
                         --avgType mean --transform obs/exp --howToCluster center --colorMap bwr &
    hicAggregateContacts --matrix results/hicexplorer/matrix/correct_all_40kb/2cell_rep1_f1_2cell_rep1_r2_40kb_nor_ICE.h5 \
                         --BED ${enhancer} --BED2 ${gene} \
                         --outFileName 2c_rep1_aggregate_Contacts_${bedpe}_0p5to1p75.pdf --vMin 0.5 --vMax 1.75 \
                         --range 1000:20000000 --numberOfBins 100 \
                         --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
                         --avgType mean --transform obs/exp --howToCluster center --colorMap bwr &
    hicAggregateContacts --matrix results/hicexplorer/matrix/correct_all_40kb/2cell_rep1_f1_2cell_rep1_r2_40kb_nor_ICE.h5 \
                         --BED ${enhancer} --BED2 ${gene} \
                         --outFileName 2c_rep1_aggregate_Contacts_${bedpe}_0p5to2p0.pdf --vMin 0.5 --vMax 2.0 \
                         --range 1000:20000000 --numberOfBins 100 \
                         --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
                         --avgType mean --transform obs/exp --howToCluster center --colorMap bwr
    hicAggregateContacts --matrix results/hicexplorer/matrix/correct_all_40kb/2cell_rep1_f1_2cell_rep1_r2_40kb_nor_ICE.h5 \
                         --BED ${enhancer} --BED2 ${gene} \
                         --outFileName 2c_rep1_aggregate_Contacts_${bedpe}_0p5to2p25.pdf --vMin 0.5 --vMax 2.25 \
                         --range 1000:20000000 --numberOfBins 100 \
                         --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
                         --avgType mean --transform obs/exp --howToCluster center --colorMap bwr &
    hicAggregateContacts --matrix results/hicexplorer/matrix/correct_all_40kb/2cell_rep1_f1_2cell_rep1_r2_40kb_nor_ICE.h5 \
                         --BED ${enhancer} --BED2 ${gene} \
                         --outFileName 2c_rep1_aggregate_Contacts_${bedpe}_0p5to2p5.pdf --vMin 0.5 --vMax 2.5 \
                         --range 1000:20000000 --numberOfBins 100 \
                         --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
                         --avgType mean --transform obs/exp --howToCluster center --colorMap bwr
    # 2 cell rep2
    hicAggregateContacts --matrix results/hicexplorer/matrix/correct_all_40kb/2cell_rep2_f1_2cell_rep2_r2_40kb_nor_ICE.h5 \
                         --BED ${enhancer} --BED2 ${gene} \
                         --outFileName 2c_rep2_aggregate_Contacts_${bedpe}_0p5to1p5.pdf --vMin 0.5 --vMax 1.5 \
                         --range 1000:20000000 --numberOfBins 100 \
                         --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
                         --avgType mean --transform obs/exp --howToCluster center --colorMap bwr &
    hicAggregateContacts --matrix results/hicexplorer/matrix/correct_all_40kb/2cell_rep2_f1_2cell_rep2_r2_40kb_nor_ICE.h5 \
                         --BED ${enhancer} --BED2 ${gene} \
                         --outFileName 2c_rep2_aggregate_Contacts_${bedpe}_0p5to1p75.pdf --vMin 0.5 --vMax 1.75 \
                         --range 1000:20000000 --numberOfBins 100 \
                         --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
                         --avgType mean --transform obs/exp --howToCluster center --colorMap bwr &
    hicAggregateContacts --matrix results/hicexplorer/matrix/correct_all_40kb/2cell_rep2_f1_2cell_rep2_r2_40kb_nor_ICE.h5 \
                         --BED ${enhancer} --BED2 ${gene} \
                         --outFileName 2c_rep2_aggregate_Contacts_${bedpe}_0p5to2p0.pdf --vMin 0.5 --vMax 2.0 \
                         --range 1000:20000000 --numberOfBins 100 \
                         --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
                         --avgType mean --transform obs/exp --howToCluster center --colorMap bwr
    hicAggregateContacts --matrix results/hicexplorer/matrix/correct_all_40kb/2cell_rep2_f1_2cell_rep2_r2_40kb_nor_ICE.h5 \
                         --BED ${enhancer} --BED2 ${gene} \
                         --outFileName 2c_rep2_aggregate_Contacts_${bedpe}_0p5to2p25.pdf --vMin 0.5 --vMax 2.25 \
                         --range 1000:20000000 --numberOfBins 100 \
                         --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
                         --avgType mean --transform obs/exp --howToCluster center --colorMap bwr &
    hicAggregateContacts --matrix results/hicexplorer/matrix/correct_all_40kb/2cell_rep2_f1_2cell_rep2_r2_40kb_nor_ICE.h5 \
                         --BED ${enhancer} --BED2 ${gene} \
                         --outFileName 2c_rep2_aggregate_Contacts_${bedpe}_0p5to2p5.pdf --vMin 0.5 --vMax 2.5 \
                         --range 1000:20000000 --numberOfBins 100 \
                         --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
                         --avgType mean --transform obs/exp --howToCluster center --colorMap bwr
    rm ${enhancer} ${gene}
done



### >>> Hic data contact frequency with ZGA genes (YHW)
# - load regions
c8_enhancer=/home/cmq/bioinfo/project-cmq/embryo_93/results/bedtools/revised_8c_enhancer/8cell_enhancer.bed
repeat_anno=/home/yhw/genome/ucsc/homo_sapiens/hg38/repeat/GRCh38_RepeatMasker_repeat_region.bed
sva_enhancer=./metadata/GRCh38_RepeatMasker_repeat_SVA_region_in_8c_enhancer.bed
bedtools intersect -a ${repeat_anno} -b ${c8_enhancer} | awk '{if(($3-$2)>=200) print $0}' > ${sva_enhancer}
k9_k27_sva=/home/cmq/bioinfo/project-cmq/embryo_93/results/idr/p_0.05/filtered/covered_repeat_age/seperate/c4_c8_93_27_shared_peaks_merged_ucsc_bl_covered_SVA_D_SVA_F_quan_50pect.bed
k9_k27_l1=/home/cmq/bioinfo/project-cmq/embryo_93/results/idr/p_0.05/filtered/covered_repeat_age/seperate/c4_c8_93_27_shared_peaks_merged_ucsc_bl_covered_L1HS_L1PA2_L1PA3_quan_50pect.bed
# - extract contact frequency from hic matrix
for file in ./tmp/*tsv
do
    prefix=`echo $file | xargs basename -s ".tsv"`
    group=`echo $file | xargs basename -s ".tsv" | sed 's,_f1,:,g' | cut -d ":" -f 1`
    interact=${file}
    for pks in sva_enhancer c8_enhancer k9_k27_sva k9_k27_l1
    do
        cat ${interact} \
           | awk -v cell=${group} -v region=${pks} 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,cell":"region}' \
           | bedtools intersect -a stdin -b ${!pks} -wa \
           | awk 'BEGIN{FS=OFS="\t"}{print $4,$5,$6,$1,$2,$3,$7,$8}' \
           | bedtools intersect -a stdin -b ./metadata/hs_zga_gene_pc_lncrna_ensembl.bed -wa > ./results/bedtools/contacts_with_zga/${group}_${pks}_contact_with_major_ZGA_gene.txt
    done
done

