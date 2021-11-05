# good luck        _      _      _
#               __(.)< __(.)> __(.)=
#               \___)  \___)  \___)
# ---------------------------------------
SHELL := /bin/bash
working_dir = $(shell pwd)
pd = /projects/marralab_cic_prj/shared_objects
h_chip = /projects/marralab_cic_prj/histone_mod_chip
chromhmm_18_states = $(shell cat ${h_chip}/chromhmm_18_states.tsv)

bedtools = /projects/sdlee_prj/installs/bedtools2/bin
bedops = /projects/jtopham_prj/installs/bedops/bin
deeptools = /gsc/software/linux-x86_64-centos6/python-2.7.12/bin
homer = /projects/sdlee_prj/installs/homer/bin
homer_jt = /projects/jtopham_prj/installs/homer/bin
finder = /home/mbilenky/bin/Solexa_Java/FindER.1.0.1e.jar
java = /gsc/software/linux-x86_64/jre1.8.0_66/bin/java
macs2 = /projects/sdlee_prj/installs/python2.7.12_venv/bin/macs2
phantom = /projects/sdlee_prj/installs/phantompeakqualtools/run_spp.R
picard_dupe = /projects/jtopham_prj/installs/picard-tools-1.114/MarkDuplicates.jar
ssp = /projects/sdlee_prj/installs/SSP/bin/ssp
UCSC = /projects/sdlee_prj/installs/UCSC_tools

lib_pairs = $(shell cat lib_pairs.tsv | awk '{print$$1"__"$$2}')
rep_pairs = $(shell cat rep_pairs.tsv | awk '{print$$1"__"$$2}')
hmod_libs = $(shell cat ${h_chip}/libs.tsv | cut -f1)

# ----------------------
#  set up directory
# ----------------------
libs ?= A77579 A77580 A77581 A77582
libs_no_input ?= A77579 A77581

setupdir = $(foreach lib, ${libs}, $(lib)/bam)
dir : ${setupdir}

# also make directories called figs/ tmp/ objects/ scripts/
# and run:
# cp /projects/vleblanc_prj/genomes/ensembl/e75_promoters_1500up_500down.bed ./tss_pad_1500_500.bed

# fix tss file to give TSS coordinates (1bp)
PERCENT := %

./tss.bed :
	cat ./tss_pad_1500_500.bed | grep "[-]1" | \
       	awk '{printf($$1);printf("${PERCENT}10.0f", ($$2 + 500)); printf("${PERCENT}10.0f", ($$3 - 1499));print("\t"$$6)}' | \
        awk '{print$$1"\t"$$2"\t"$$3"\t"$$4}' > ./tmp/tmp.bed; \
       	cat ./tss_pad_1500_500.bed | grep -v "[-]1" | \
        awk '{printf($$1);printf("${PERCENT}10.0f", ($$2 + 1500)); printf("${PERCENT}10.0f", ($$3 - 499));print("\t"$$6)}' | \
	awk '{print$$1"\t"$$2"\t"$$3"\t"$$4}' | awk '$$2 > 0' | cat ./tmp/tmp.bed - | \
       	sort-bed --max-mem 10G --tmpdir ./tmp - > $@

# move on to processing bams
%/bam :
	mkdir $*; mkdir $*/bam;	mkdir $*/peaks;

# ----------------------
#  process BAM files
# ----------------------
# note this was all done using samtools 0.1.18

flagstat_loop = $(foreach lib, ${libs}, $(lib)/bam/flagstat.txt)
flagstat : ${flagstat_loop}

%/bam/flagstat.txt :
	samtools flagstat $*/bam/gsc.bam > $@

# filter unmapped reads; sort and index
filter_sort_loop = $(foreach lib, ${libs}, $(lib)/bam/gsc.sorted.bam)
filter_sort : ${filter_sort_loop}

%/bam/gsc.sorted.bam :
	samtools view -bF 4 $*/bam/gsc.bam > $*/bam/tmp.bam; \
	samtools sort $*/bam/tmp.bam $*/bam/gsc.sorted; \
	samtools index $*/bam/gsc.sorted.bam $*/bam/gsc.sorted.bam.bai; rm $*/bam/tmp.bam

rmdup = $(foreach lib, ${libs}, $(lib)/bam/gsc.nodupe.bam)
rmdupe : ${rmdup}

# make files with dupes removed, for findER
%/bam/gsc.nodupe.bam :
	java -jar ${picard_dupe} I=$*/bam/gsc.sorted.bam O=$*/bam/tmp.bam M=$*/bam/dupe.stdout.txt \
	ASSUME_SORTED=true REMOVE_DUPLICATES=true; \
	samtools sort $*/bam/tmp.bam $*/bam/gsc.nodupe; \
	samtools index $@ $@.bai; rm $*/bam/tmp.bam

# filter multi-mapping reads (https://wabi-wiki.scilifelab.se/display/KB/Filter+uniquely+mapped+reads+from+a+BAM+file)
filter_unique_loop = $(foreach lib, ${libs}, $(lib)/bam/gsc.sorted.unique.bam)
filter_unique : ${filter_unique_loop}

%/bam/gsc.sorted.unique.bam :
	samtools view -b -q 5 -F 256 $*/bam/gsc.sorted.bam > $*/bam/tmp.bam; \
	samtools sort $*/bam/tmp.bam $*/bam/gsc.sorted.unique; \
	samtools index $*/bam/gsc.sorted.unique.bam $*/bam/gsc.sorted.unique.bam.bai; \
	rm $*/bam/tmp.bam

A77582/bam/gsc.sorted.unique.bam :
	samtools view -b -q 5 -F 256 A77582/bam/gsc.sorted.bam > A77582/bam/tmp.bam; \
	samtools sort A77582/bam/tmp.bam A77582/bam/gsc.sorted.unique; \
	samtools index A77582/bam/gsc.sorted.unique.bam A77582/bam/gsc.sorted.unique.bam.bai

# ------------------------------------------------------------------------------
#  get library fragment sizes
# ------------------------------------------------------------------------------
frag_size_loop = $(foreach lib, ${libs}, ./$(lib)/bam/fragment_sizes.txt)
frag_sizes : ${frag_size_loop}

%/bam/fragment_sizes.txt :
	samtools view -f 0x02 $*/bam/gsc.bam | awk '{if ($$9 > 0) {print $$9} }' > $@

# ------------------------------------------------------------------------------
#  Generate BigWig Files (bin size = 10, normalize to RPKM)
# ------------------------------------------------------------------------------
bigwig_loop = $(foreach lib, ${libs}, $(lib)/bam/gsc.nodupe.rpkm.bw)
bigwig : ${bigwig_loop}

%/bam/gsc.nodupe.rpkm.bw :
	${deeptools}/bamCoverage -p 70 -b $*/bam/gsc.nodupe.bam -o $@ \
	--binSize 10 \
	--normalizeUsing RPKM \
	--ignoreForNormalization chrX \
	--extendReads

# subtract input rpkm from ip rpkm
bw_subtract_loop = $(foreach pair, ${lib_pairs}, $(pair)/bam/gsc.nodupe.subtract.bw)
bw_subtract : ${bw_subtract_loop}

%/bam/gsc.nodupe.subtract.bw :
	ip=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	input=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	${deeptools}/bigwigCompare -p 60 \
	-b1 $${ip}/bam/gsc.nodupe.rpkm.bw \
	-b2 $${input}/bam/gsc.nodupe.rpkm.bw \
	--binSize 10 \
	--operation subtract \
	--skipNAs \
	-o $${ip}/bam/gsc.nodupe.subtract.bw

# normalize to input
bw_log2ratio_loop = $(foreach pair, ${lib_pairs}, $(pair)/bam/gsc.sorted.log2ratio.bw)
bw_log2ratio : ${bw_log2ratio_loop}

%/bam/gsc.sorted.log2ratio.bw :
	ip=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	input=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	${deeptools}/bamCompare -p 20 \
	-b1 $${ip}/bam/gsc.sorted.bam \
	-b2 $${input}/bam/gsc.sorted.bam \
	--scaleFactorsMethod SES \
	--ignoreForNormalization chrX \
	--extendReads \
	-o $${ip}/bam/gsc.sorted.log2ratio.bw

ucsc_loop = $(foreach lib, ${libs}, $(lib)/ucsc_bam.bw)
ucsc : ${ucsc_loop}

%/ucsc_bam.bw :
	cp $*/bam/gsc.sorted.bam.bw /gsc/www/bcgsc.ca/downloads/sdlee/CIC/first_chip_sub/$*.bam.sorted.bw; \
	echo "track type=bigWig name='$*' description='$*' \
	color=70,130,180 bigDataUrl=http://bigwigviewer:bwvbwv@www.bcgsc.ca/downloads/sdlee/CIC/first_chip_sub/$*.bam.sorted.bw" >> upload

# ------------------------------------------------------------------------------
#  GC Bias
# ------------------------------------------------------------------------------
check_GC = $(foreach lib, ${libs}, $(lib)/bam/GCfreq.txt)
check_GC_bias : ${check_GC}

# generated by using faToTwoBit (downloaded from http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/)
# on human_g1k_v37.fasta (in /projects/marralab_cic_prj/shared_objects/reference_genome/)
2bit = /projects/marralab_cic_prj/shared_objects/reference_genome/human_g1k_v37/human_g1k_v37.2bit
#	--effectiveGenomeSize 2747877777 \

%/bam/GCfreq.txt :
	${deeptools}/computeGCBias -p 20 -b $*/bam/gsc.sorted.unique.bam \
	--effectiveGenomeSize 2747877777 \
	-g ${2bit} \
	-o $@ \
	--biasPlot $*/bam/GCbiasplot.png

GC_pairs = $(shell cat GC_pairs.tsv | awk '{print$$1"__"$$2}')
correct_GC_loop = $(foreach pair, ${GC_pairs}, $(pair)/gsc.sorted.gccorr.bam)
correct_GC : ${correct_GC_loop}

%/gsc.sorted.gccorr.bam :
	bam=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	input=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	${deeptools}/correctGCBias -p 20 -b $${bam}/bam/gsc.sorted.unique.bam \
	--effectiveGenomeSize 2747877777 \
	-g ${2bit} \
	-freq $${input}/bam/GCfreq.txt \
	-o $${bam}/bam/gsc.sorted.gccorr.bam

# ------------------------------------------------------------------------------
#  Correlate Bams and Plot
# ------------------------------------------------------------------------------
bams = $(shell ls ./*/bam/gsc.bam)
gccorr_bams = $(shell ls ./*/bam/gsc.sorted.gccorr.bam)

corr_bams :
	${deeptools}/multiBamSummary bins \
	--bamfiles ${bams} \
	--minMappingQuality 30 \
	--labels ${libs} \
	-out objects/readCounts.npz --outRawCounts objects/readCounts.tab

corr_spearman_heatmap :
	${deeptools}/plotCorrelation \
	-in objects/readCounts.npz \
	--corMethod spearman --skipZeros \
	--plotTitle "Spearman Correlation of Read Counts" \
	--whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
	-o figs/QC_SL/heatmap_spearman_readCounts.png \
	--outFileCorMatrix objects/spearman_readCounts.tab

correlate_bams : corr_bams corr_spearman_heatmap

# ------------------------------------------------------------------------------
#  Correlate BigWigs and plot
# ------------------------------------------------------------------------------
# Preferred over MultiBamSummary since you compare normalized read coverages

bigwigs = $(shell ls ./*/bam/gsc.sorted.rpkm.bw)

#bins are 10kb by default
corr_bw :
	${deeptools}/multiBigwigSummary bins \
	-p 10 \
	--bwfiles ${bigwigs} \
	--labels NHA_rep1_chip NHA_rep1_input NHA_rep2_chip NHA_rep2_input \
	-out objects/bw_readCounts.npz \
	--outRawCounts objects/bw_readCounts.tab

corr_bw_sp_heatmap :
	${deeptools}/plotCorrelation \
	-in objects/bw_readCounts.npz \
	--corMethod spearman --skipZeros \
	--plotTitle "Spearman Correlation of Normalized Read Counts (RPKM)" \
	--whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
	-o figs/QC_SL/heatmap_sp_readCounts.png \
	--outFileCorMatrix objects/spearman_readCounts.tab

correlate_bw : corr_bw corr_bw_sp_heatmap

# ------------------------------------------------------------------------------
#  Principal Component Analysis
# ------------------------------------------------------------------------------
PCA :
	${deeptools}/plotPCA \
	-in objects/bw_readCounts.npz \
	-o figs/QC_SL/PCA_bw_readCounts.png \
	-T "PCA of Normalized Read Counts" \
	--outFileNameData objects/PCA_bw.tab

# ------------------------------------------------------------------------------
#  Fingerprint
# ------------------------------------------------------------------------------
fingerprint :
	${deeptools}/plotFingerprint \
	-b ${bams} -p10 \
	--extendReads \
	--labels NHA_rep1_chip NHA_rep1_input NHA_rep2_chip NHA_rep2_input \
	--minMappingQuality 30 --skipZeros \
	-T "Fingerprints of First Submission Samples"  \
	--plotFile figs/QC_SL/fingerprints.png \
	--outRawCounts objects/fingerprints.tab

fingerprint_gccorr :
	${deeptools}/plotFingerprint \
	-b ${gccorr_bams} -p10 \
	--extendReads \
	--labels NHA_rep1_chip NHA_rep1_input NHA_rep2_chip NHA_rep2_input \
	--minMappingQuality 30 --skipZeros \
	-T "Fingerprints of First Submission Samples"  \
	--plotFile figs/QC_SL/fingerprints.png \
	--outRawCounts objects/fingerprints.tab

# ------------------------------------------------------------------------------
#  Cross Correlation Analysis
# ------------------------------------------------------------------------------
# added library(caTools) to ${phantom}, otherwise get an error "cannot find function runmean"

phantom_loop = $(foreach lib, ${libs}, $(lib)/bam/phantom.txt)
phantom : ${phantom_loop}

%/bam/phantom.txt :
	Rscript ${phantom} \
	-p=20 \
	-s=0:5:500 \
	-c=$*/bam/gsc.sorted.unique.bam \
	-savp \
	-out=$@

ssp_loop = $(foreach lib, ${libs}, $(lib)/bam/ssp)
ssp : ${ssp_loop}

%/bam/ssp :
	${ssp} -p10 \
	-i $*/bam/gsc.sorted.bam \
	-o ssp \
	--odir $*/bam/ssp \
	--pair \
	--gt ${pd}hg19_chromsizes.txt

# ------------------------------------------------------------------------------
#  call peaks using findER
# ------------------------------------------------------------------------------
finder_loop = $(foreach lib_pair, ${lib_pairs}, $(lib_pair)/peaks.bed)
finder : ${finder_loop}

%/peaks.bed :
	ip=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	input=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	${java} -jar -Xmx20G ${finder} -maxGap 0 \
	-signalBam $${ip}/bam/gsc.nodupe.bam \
        -inputBam $${input}/bam/gsc.nodupe.bam \
	-out $${ip}/peaks/; \
        gunzip -c $${ip}/peaks/gsc.nodupe.vs.gsc.nodupe.FDR_0.05.FindER.bed.gz | \
        awk '{gsub("chr", "");print}' | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - \
	> ./$${ip}/peaks/peaks.sorted.bed; \
        rm ./$${ip}/peaks/gsc.nodupe.vs.gsc.nodupe.FDR_0.05.FindER.bed.gz

# ------------------------------------------------------------------------------
#  call peaks using MACS2
# ------------------------------------------------------------------------------
macs2_loop = $(foreach lib_pair, ${lib_pairs}, $(lib_pair)/macs2.bed)
macs2 : ${macs2_loop}

# effective genome size for hg19 (75bp readlen) = 2736124973
%/macs2.bed :
	ip=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	input=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	${macs2} callpeak -B \
	-t $${ip}/bam/gsc.sorted.bam \
	-c $${input}/bam/gsc.sorted.bam \
	-f BAMPE \
	-g 2736124973 \
	-q 0.01 \
	-n macs2_q0.01 \
	--outdir $${ip}/peaks

macs2_0.05_loop = $(foreach lib_pair, ${lib_pairs}, $(lib_pair)/macs2.bed)
macs2_0.05 : ${macs2_0.05_loop}

# effective genome size for hg19 (75bp readlen) = 2736124973
%/macs2.bed :
	ip=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	input=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	${macs2} callpeak -B \
	-t $${ip}/bam/gsc.sorted.bam \
	-c $${input}/bam/gsc.sorted.bam \
	-f BAMPE \
	-g 2736124973 \
	-q 0.05 \
	-n macs2_q0.05 \
	--outdir $${ip}/peaks

# ------------------------------------------------------------------------------
#  filter blacklisted regions from CIC peaks
# ------------------------------------------------------------------------------
# cat ../shared_objects/consensusBlacklist.bed | sed 's/^chr//' > tmp/bl.bed
A77581/peaks/mq0.05.rmbl.bed :
	${bedtools}/subtractBed -A -a A77581/peaks/macs2_q0.05_peaks.narrowPeak \
	-b ../histone_mod_chip/tmp/bl.bed > $@

# ------------------------------------------------------------------------------
#  Obtain RPKM within MACS2 peaks
# ------------------------------------------------------------------------------
mq0.05_rpkm_loop = $(foreach pair, ${lib_pairs}, $(pair)/peaks/mq0.05_rpkm.bed)
mq0.05_rpkm : ${mq0.05_rpkm_loop}

%/peaks/mq0.05_rpkm.bed :
	ip=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	input=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	${deeptools}/multiBigwigSummary BED-file -p 60 \
	-b $${ip}/bam/gsc.nodupe.rpkm.bw $${input}/bam/gsc.nodupe.rpkm.bw \
	--BED $${ip}/peaks/macs2_q0.05_peaks.narrowPeak \
	--labels chip input \
	--outRawCounts $${ip}/peaks/mq0.05_rpkm.tab \
	-o $${ip}/peaks/mq0.05_rpkm.npz; \
	awk 'NR>1' $${ip}/peaks/mq0.05_rpkm.tab | sort -k1,1 -k2,2n > $${ip}/peaks/mq0.05_rpkm.bed

# ------------------------------------------------------------------------------
#  annotate MACS2 peaks to nearest tss
# ------------------------------------------------------------------------------
# Annotate peaks with nearest tss of protein coding genes
anno_peak_pc_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/mq0.05_tss_pc.bed)
anno_peak_pc : ${anno_peak_pc_loop}

%/peaks/mq0.05_tss_pc.bed :
	${bedtools}/closestBed -D "ref" -a $*/peaks/macs2_q0.05_peaks.narrowPeak -b ${pd}/tss.pc.bed > $@


# ------------------------------------------------------------------------------
#  coverage in MACS2 peaks
# ------------------------------------------------------------------------------
./A77581/peaks/mq0.01_cov.npz :
	${deeptools}/multiBigwigSummary BED-file -p 40 \
	-b ./A77581/bam/gsc.sorted.rpkm.bw ./A77582/bam/gsc.sorted.rpkm.bw \
	--BED ./A77581/peaks/macs2_q0.01_peaks.narrowPeak \
	--labels NHA_rep2_chip NHA_rep2_input \
	--outRawCounts ./A77581/peaks/mq0.01_cov.tab \
	-o $@

# ------------------------------------------------------------------------------
#  Motif Analysis
# ------------------------------------------------------------------------------
# motif analysis on top 150 filtered peaks, centred on the summit
/motif/homer/coenriched_summit_hc/knownResults.html :
	cat NHA_CIC_hc_peaks.bed | \
	/projects/jtopham_prj/installs/homer/bin/findMotifsGenome.pl - \
	${pd}/hg19.fa motifs/homer/coenriched_summit_hc/ -p 40 -S 20 -size 200 -mask

# high confidence CIC peak set by merging rep2 consensus peaks with peaks found in G144 from Weissmann et al
# ------------------------------------------------------------------------------
#  intersect Weissman et al G144
# ------------------------------------------------------------------------------
./rep2_i_G144/intersect_peaks.bed :
	cut -f1-3 A77581/peaks/macs2_q0.01_peaks.narrowPeak | \
	awk '$$1 !~ /GL/' > A77581/peaks/mq0.01.bed; \
	${bedtools2}/intersectBed -wo -a A77581/peaks/mq0.01.bed \
	-b ../Weissman_et_al/paper_peaks/G144_24h_MEKi.hg19.sorted.bed > $@

# generate matrix for CIC peaks (mq0.01 rep2)
./A77581/peaks/mq0.01_peak_mat.gz :
	${deeptools}/computeMatrix reference-point -p 40 \
	--referencePoint center \
	-b 5000 -a 5000 \
	-R ./A77581/peaks/peaks.sorted.bed \
	-S ./A77581/bam/gsc.sorted.rpkm.bw \
	-o $@

# pad mq0.01 peaks by 5kb from peak center and bin by 100bp
./A77581/peaks/mq0.01_5kb.bed :
	cat ./A77581/peaks/mq0.01.bed | \
	awk '{print $$1"\t"int(($$2 + $$3)/2)}' | \
	awk '{print $$1"\t"($$2 - 5000)"\t"($$2 + 5000)}' | \
	awk '$$2 > 0' | \
	${bedtools2}/windowMaker -b - -w 100 -i srcwinnum | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $@

marks ?= h3k4me1 h3k4me3 h3k27ac h3k27me3 h3k36me3 h3k9me3

# obtain rpkm for CIC and histone mod ChIP across padbinned mq0.01 peaks
mq0.01_rep2_hcov_loop = $(foreach mark, ${marks}, ./A77581/peaks/mq0.01_5kb_$(mark)_rpkm.tab)
mq0.01_rep2_hcov : ${mq0.01_rep2_hcov_loop}

./A77581/peaks/mq0.01_5kb_%_rpkm.tab :
	${deeptools}/multiBigwigSummary BED-file -p 60 \
	--bwfiles ./A77581/bam/gsc.sorted.rpkm.bw ../histone_mod_chip/merged_bw/*$*.rpkm.bw \
	--BED ./A77581/peaks/mq0.01_5kb.bed \
	--outRawCounts $@ \
	-out ./A77581/peaks/mq0.01_5kb_$*_rpkm.npz

rep2_i_G144_hcov_loop = $(foreach mark, ${marks}, ./rep2_i_G144/$(mark)/union_mat.gz)
rep2_i_G144_hcov : ${rep2_i_G144_hcov_loop}

./rep2_i_G144/%/union_mat.gz :
	mkdir rep2_i_G144/$*; \
	samples=$$(echo ../histone_mod_chip/*$*/bam/gsc.sorted.rpkm.bw | awk '{ gsub("../histone_mod_chip/", ""); print }' | \
	awk '{ gsub("/bam/gsc.sorted.rpkm.bw", ""); print }'); \
	${deeptools}/computeMatrix reference-point -p 20 \
	--referencePoint center \
	-b 10000 -a 10000 \
	-R rep2_i_G144/union_peaks.bed \
	-S ../histone_mod_chip/*$*/bam/gsc.sorted.rpkm.bw \
	-o $@

# final peak set obtained in R (overlap between both replicates and/or between G144 from Weissmann et al)
# ------------------------------------------------------------------------------
#  coverage within final peak set
# ------------------------------------------------------------------------------
final_peakcov_loop = $(foreach lib_pair, ${lib_pairs}, $(lib_pair)/cov_in_final_peaks.bed)
final_peakcov : ${final_peakcov_loop}

%/cov_in_final_peaks.bed :
	ip=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	input=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	${bedtools}/coverageBed -abam $${ip}/bam/gsc.sorted.bam -b \
	$${ip}/peaks/final_peaks.bed -counts | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $${ip}/peaks/cov_in_final_peaks_tmp.bed; \
	${bedtools}/coverageBed -abam $${input}/bam/gsc.sorted.bam -b \
	$${ip}/peaks/peaks.sorted.bed -counts | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $${input}/peaks/cov_in_final_peaks.bed; \
	paste $${ip}/peaks/cov_in_final_peaks_tmp.bed <(awk '{print $$5}' $${input}/peaks/cov_in_final_peaks.bed) > \
	$${ip}/peaks/cov_in_final_peaks.bed

# ------------------------------------------------------------------------------
#  histone coverage within final peak set
# ------------------------------------------------------------------------------
final_peaks_hcov_loop = $(foreach lib, ${hmod_libs}, ./final_peakcov/$(lib)_cov_final_peaks.bed)
final_peaks_hcov : ${final_peaks_hcov_loop}

./final_peakcov/%_cov_final_peaks.bed :
	${bedtools}/coverageBed -abam ${h_chip}/$*/bam/gsc.sorted.bam \
	-b ./merged_peaks/final_peaks_2.bed -counts | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $@

# ------------------------------------------------------------------------------
#  center peaks around MACS2 summits, fix width to 500bp and obtain sequence
# ------------------------------------------------------------------------------
summit_500_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/cons_summit_500.bed)
summit_500 : ${summit_500_loop}

%/peaks/cons_summit_500.bed :
	cat $*/peaks/cons_peaks.bed | \
	awk '{print $$1"\t"($$2 + $$10)}' | \
	awk '{print $$1"\t"($$2 - 250)"\t"($$2 + 250)}' | \
	sort -u > $@; \
	${bedtools}/fastaFromBed -fi ${pd}/hg19.fa -bed $@ -fo $*/peaks/cons_summit_500.fa

summit_both_500_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/cons_both_summit_500.bed)
summit_both_500 : ${summit_both_500_loop}

%/peaks/cons_both_summit_500.bed :
	cat $*/peaks/cons_peaks_both.bed | \
	awk '{print $$1"\t"($$2 + $$5)}' | \
	awk '{print $$1"\t"($$2 - 250)"\t"($$2 + 250)}' | \
	sort -u > $@; \
	${bedtools}/fastaFromBed -fi ${pd}/hg19.fa -bed $@ -fo $*/peaks/cons_both_summit_500.fa

final_summit_500_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/final_summit_500.bed)
final_summit_500 : ${final_summit_500_loop}

%/peaks/final_summit_500.bed :
	cat $*/peaks/final_peaks.bed | \
	awk '{print $$1"\t"($$2 + $$5)}' | \
	awk '{print $$1"\t"($$2 - 250)"\t"($$2 + 250)}' | \
	sort -u > $@; \
	${bedtools}/fastaFromBed -fi ${pd}/hg19.fa -bed $@ -fo $*/peaks/final_summit_500.fa

final_summit_100_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/final_summit_100.bed)
final_summit_100 : ${final_summit_100_loop}

%/peaks/final_summit_100.bed :
	cat $*/peaks/final_peaks.bed | \
	awk '{print $$1"\t"($$2 + $$5)}' | \
	awk '{print $$1"\t"($$2 - 50)"\t"($$2 + 50)}' | \
	sort -u > $@; \
	${bedtools}/fastaFromBed -fi ${pd}/hg19.fa -bed $@ -fo $*/peaks/final_summit_100.fa

# ------------------------------------------------------------------------------
#  fix width around peak center and obtain sequence
# ------------------------------------------------------------------------------
final_center_500_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/final_center_500.bed)
final_center_500 : ${final_center_500_loop}

%/peaks/final_center_500.bed :
	cat $*/peaks/final_peaks.bed | \
	awk '{print $$1"\t"(int(($$2 + $$3) / 2) + 0.5)}' | \
	awk '{print $$1"\t"($$2 - 250)"\t"($$2 + 250)}' | \
	sort -u > $@; \
	${bedtools}/fastaFromBed -fi ${pd}/hg19.fa -bed $@ -fo $*/peaks/final_center_500.fa

final_center_100_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/final_center_100.bed)
final_center_100 : ${final_center_100_loop}

%/peaks/final_center_100.bed :
	cat $*/peaks/final_peaks.bed | \
	awk '{print $$1"\t"(int(($$2 + $$3) / 2) + 0.5)}' | \
	awk '{print $$1"\t"($$2 - 250)"\t"($$2 + 250)}' | \
	sort -u > $@; \
	${bedtools}/fastaFromBed -fi ${pd}/hg19.fa -bed $@ -fo $*/peaks/final_center_100.fa

# ------------------------------------------------------------------------------
#  Identify enriched motifs in CIC consensus peaks
# ------------------------------------------------------------------------------
cons_peaks_motifs_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/knownResults.html)
cons_peaks_motifs : ${cons_peaks_motifs_loop}

%/peaks/knownResults.html :
	cat $*/peaks/cons_peaks.bed | \
	/projects/jtopham_prj/installs/homer/bin/findMotifsGenome.pl - \
	/projects/jtopham_prj/Thesis_prj/data/HEK_RNAseq/hg19.fa $*/peaks/homer_100bp -p 20 -size 50

#in rep1 & 2 peaks that overlap
cons_peaks_both_motifs_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/both_homer_100bp/knownResults.html)
cons_peaks_both_motifs : ${cons_peaks_both_motifs_loop}

%/peaks/both_homer_100bp/knownResults.html :
	cat $*/peaks/cons_peaks_both.bed | \
	/projects/jtopham_prj/installs/homer/bin/findMotifsGenome.pl - \
	/projects/jtopham_prj/Thesis_prj/data/HEK_RNAseq/hg19.fa $*/peaks/both_homer_100bp -p 20 -size 50

final_peaks_motifs_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/final_homer_100bp/knownResults.html)
final_peaks_motifs : ${final_peaks_motifs_loop}

%/peaks/final_homer_100bp/knownResults.html :
	cat $*/peaks/final_summit_500.bed | \
	/projects/jtopham_prj/installs/homer/bin/findMotifsGenome.pl - \
	/projects/jtopham_prj/Thesis_prj/data/HEK_RNAseq/hg19.fa $*/peaks/final_homer_100bp -p 20 -size 50

final_center_motifs_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/final_center_homer/knownResults.html)
final_center_motifs : ${final_center_motifs_loop}

%/peaks/final_center_homer/knownResults.html :
	cat $*/peaks/final_center_500.bed | \
	/projects/jtopham_prj/installs/homer/bin/findMotifsGenome.pl - \
	/projects/jtopham_prj/Thesis_prj/data/HEK_RNAseq/hg19.fa $*/peaks/final_center_homer -p 20 -size 50

# ------------------------------------------------------------------------------
#  Padbin consensus peaks centred on MACS2 summit
# ------------------------------------------------------------------------------
padbin_cons_summit_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/cons_peaks_padbin_2kb.bed)
padbin_cons_summit : ${padbin_cons_summit_loop}

%/peaks/cons_peaks_padbin_2kb.bed :
	cat $*/peaks/cons_peaks.bed | \
	awk '{print $$1"\t"($$2 + $$10)}' | \
	awk '{print $$1"\t"($$2 - 2000)"\t"($$2 + 2000)}' | \
	awk '$$2 > 0' | \
	${bedtools}/windowMaker -b - -w 100 -i srcwinnum | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $@

# ------------------------------------------------------------------------------
#  Padbin final peaks centred on MACS2 peak center
# ------------------------------------------------------------------------------
padbin_final_peaks_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/final_peaks_padbin.bed)
padbin_final_peaks : ${padbin_final_peaks_loop}

%/peaks/final_peaks_padbin.bed :
	cat $*/peaks/final_peaks.bed | \
	awk '{print $$1"\t"(int(($$2 + $$3)/2) + 0.5)}' | \
	awk '{print $$1"\t"($$2 - 2000)"\t"($$2 + 2000)}' | \
	awk '$$2 > 0' | \
	${bedtools}/windowMaker -b - -w 100 -i srcwinnum | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $@

cov_padbin_final_peak_center_loop = $(foreach lib_pair, ${lib_pairs}, $(lib_pair)/cov_padbin_finder_peak_center.bed)
cov_padbin_final_peak_center : ${cov_padbin_final_peak_center_loop}

%/cov_padbin_finder_peak_center.bed :
	ip=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	input=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	${bedtools}/coverageBed -abam $${ip}/bam/gsc.sorted.bam \
	-b $${ip}/peaks/final_peaks_padbin.bed -counts | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $${ip}/peaks/final_peaks_padbin_cov_tmp.bed; \
	${bedtools}/coverageBed -abam $${input}/bam/gsc.sorted.bam \
	-b $${ip}/peaks/final_peaks_padbin.bed -counts | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $${input}/peaks/final_peaks_padbin_cov.bed; \
	paste $${ip}/peaks/final_peaks_padbin_cov_tmp.bed <(awk '{print $$5}' $${input}/peaks/final_peaks_padbin_cov.bed) > \
	$${ip}/peaks/final_peaks_padbin_cov.bed

# ------------------------------------------------------------------------------
#  overlap peaks between Weissman et al 2018
# ------------------------------------------------------------------------------
intersect_Weissman_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/mq0.01_i_hNSC18.5.bed)
intersect_Weissman : ${intersect_Weissman_loop}

%/peaks/mq0.01_i_hNSC18.5.bed :
	${bedtools}/intersectBed -wo -a $*/peaks/macs2_q0.01_peaks.narrowPeak \
	-b ../Weissman_et_al/paper_peaks/hNSC18.5_24h_MEKi.hg19.sorted.bed > \
	$*/peaks/mq0.01_i_hNSC18.5.bed; \
	${bedops}/closest-features --dist --closest --delim '\t' \
	$*/peaks/mq0.01_i_hNSC18.5.bed ${pd}/tss_pad_1500_500.bed > \
	$*/peaks/mq0.01_i_hNSC18.5_tss.bed; \

	${bedtools}/intersectBed -wo -a $*/peaks/macs2_q0.01_peaks.narrowPeak \
	-b ../Weissman_et_al/paper_peaks/hNSC21.5_24h_MEKi.hg19.sorted.bed > \
	$*/peaks/mq0.01_i_hNSC21.5.bed; \
	${bedops}/closest-features --dist --closest --delim '\t' \
	$*/peaks/mq0.01_i_hNSC21.5.bed ${pd}/tss_pad_1500_500.bed > \
	$*/peaks/mq0.01_i_hNSC21.5_tss.bed; \

	${bedtools}/intersectBed -wo -a $*/peaks/macs2_q0.01_peaks.narrowPeak \
	-b ../Weissman_et_al/paper_peaks/G144_24h_MEKi.hg19.sorted.bed > \
	$*/peaks/mq0.01_i_G144.bed; \
	${bedops}/closest-features --dist --closest --delim '\t' \
	$*/peaks/mq0.01_i_G144.bed ${pd}/tss_pad_1500_500.bed > \
	$*/peaks/mq0.01_i_G144_tss.bed

f_intersect_Weissman_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/f_i_hNSC18.5.bed)
f_intersect_Weissman : ${f_intersect_Weissman_loop}

%/peaks/f_i_hNSC18.5.bed :
	${bedtools}/intersectBed -wo -a $*/peaks/peaks.sorted.bed \
	-b ../Weissman_et_al/paper_peaks/hNSC18.5_24h_MEKi.hg19.sorted.bed > \
	$*/peaks/f_i_hNSC18.5.bed; \
	${bedops}/closest-features --dist --closest --delim '\t' \
	$*/peaks/f_i_hNSC18.5.bed ${pd}/tss_pad_1500_500.bed > \
	$*/peaks/f_i_hNSC18.5_tss.bed; \

	${bedtools}/intersectBed -wo -a $*/peaks/peaks.sorted.bed \
	-b ../Weissman_et_al/paper_peaks/hNSC21.5_24h_MEKi.hg19.sorted.bed > \
	$*/peaks/f_i_hNSC21.5.bed; \
	${bedops}/closest-features --dist --closest --delim '\t' \
	$*/peaks/f_i_hNSC21.5.bed ${pd}/tss_pad_1500_500.bed > \
	$*/peaks/f_i_hNSC21.5_tss.bed; \

	${bedtools}/intersectBed -wo -a $*/peaks/peaks.sorted.bed \
	-b ../Weissman_et_al/paper_peaks/G144_24h_MEKi.hg19.sorted.bed > \
	$*/peaks/f_i_G144.bed; \
	${bedops}/closest-features --dist --closest --delim '\t' \
	$*/peaks/f_i_G144.bed ${pd}/tss_pad_1500_500.bed > \
	$*/peaks/f_i_G144_tss.bed

#intersect consensus peaks with Weissman
cons_i_Weissman_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/cons_i_G144.bed)
cons_i_Weissman : ${cons_i_Weissman_loop}

%/peaks/cons_i_G144.bed :
	${bedtools}/intersectBed -wo -a $*/peaks/cons_peaks.bed \
	-b ../Weissman_et_al/paper_peaks/G144_24h_MEKi.hg19.sorted.bed | \
	awk '{print $$1"\t"$$2"\t"$$3"\t"$$4"\t"$$10"\t"$$6}' | sort -u > $@; \
	${bedtools}/intersectBed -wo -a $*/peaks/cons_peaks_both.bed \
	-b ../Weissman_et_al/paper_peaks/G144_24h_MEKi.hg19.sorted.bed | \
	awk '{print $$1"\t"$$2"\t"$$3"\t"$$4"\t"$$10"\t"$$6}' | sort -u > $*/peaks/cons_both_i_G144.bed

merge_final_peaks :
	cat ./merged_peaks/final_peaks.bed | \
	sort -k1,1 -k2,2n | \
	sed 's/^chr//' | \
	${bedtools}/mergeBed > ./merged_peaks/final_peaks_merged.bed

# ------------------------------------------------------------------------------
#  Convert RPGC normalized BigWig into BedGraph
# ------------------------------------------------------------------------------
bw_to_bg_loop = $(foreach lib, ${libs}, $(lib)/bam/gsc.sorted.rpgc.bedGraph)
bw_to_bg : ${bw_to_bg_loop}

%/bam/gsc.sorted.rpgc.bedGraph :
	${UCSC}/bigWigToBedGraph $*/bam/gsc.sorted.bam.bw $@

macs2_loop = $(foreach lib_pair, ${lib_pairs}, $(lib_pair)/macs2.bed)
macs2 : ${macs2_loop}

# ------------------------------------------------------------------------------
#  Call peaks on RPGC normalized BedGraph using MACS2
# ------------------------------------------------------------------------------
mq0.01_rpgc_loop = $(foreach lib_pair, ${lib_pairs}, $(lib_pair)/mq0.01_rpgc.bed)
mq0.01_rpgc : ${mq0.01_rpgc_loop}

%/mq0.01_rpgc.bed :
	ip=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	input=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	${macs2} bdgpeakcall -B \
	-t $${ip}/bam/gsc.sorted.rpgc.bedGraph \
	-c $${input}/bam/gsc.sorted.rpgc.bedGraph \
	-g 2736124973 \
	-q 0.01 \
	-n mq0.01_rpgc \
	--outdir $${ip}/peaks

# ------------------------------------------------------------------------------
#  Obtain histone mod chip coverage around CIC (f_i_G144) peaks
# ------------------------------------------------------------------------------
./merged_peaks/control_padbin_2kb.bed :
	cat merged_peaks/rep2_f_i_G144_c.bed | \
	awk '{print $$1"\t"(($$2 + $$3)/2)}' | \
	awk '{print $$1"\t"($$2 - 2000)"\t"($$2 + 2000)}' | \
	awk '$$2 > 0' | \
	${bedtools}/windowMaker -b - -w 100 -i srcwinnum | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $@

hmod_cov_loop = $(foreach lib, ${hmod_libs}, ./merged_peaks/f_i_G144_padbin_$(lib)_cov.bed)
hmod_cov : ${hmod_cov_loop}

./merged_peaks/f_i_G144_padbin_%_cov.bed :
	${bedtools}/coverageBed -abam ${h_chip}/$*/bam/gsc.sorted.bam \
	-b merged_peaks/rep2_f_i_G144_padbin_2kb.bed -counts | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $@

hmod_cov_cons_peaks_loop = $(foreach lib, ${hmod_libs}, ./merged_peaks/cons_peaks_$(lib)_cov.bed)
hmod_cov_cons_peaks : ${hmod_cov_cons_peaks_loop}

./merged_peaks/cons_peaks_%_cov.bed :
	${bedtools}/coverageBed -abam ${h_chip}/$*/bam/gsc.sorted.bam \
	-b A77579/peaks/cons_peaks_padbin_2kb.bed -counts | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > ./merged_peaks/cons_peaks_$*_cov_rep1.bed; \
	${bedtools}/coverageBed -abam ${h_chip}/$*/bam/gsc.sorted.bam \
	-b A77581/peaks/cons_peaks_padbin_2kb.bed -counts | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > ./merged_peaks/cons_peaks_$*_cov_rep2.bed

# ------------------------------------------------------------------------------
#  Generate a set of control peaks wth same size distribution as f_i_G144
# ------------------------------------------------------------------------------
merged_peaks/rep2_f_i_G144_c.bed :
	${bedtools2}/shuffleBed -i merged_peaks/rep2_f_i_G144.bed \
	-g ${pd}/hg19_chromsizes.txt \
	-excl ${pd}/consensusBlacklist.bed \
	-noOverlapping \
	-seed 12072013 \
	> $@

# control peaks final 500bp
make_ctrl_peaks_500_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/final_summit_500_ctrl.bed)
make_ctrl_peaks_500 : ${make_ctrl_peaks_500_loop}

%/peaks/final_summit_500_ctrl.bed :
	${bedtools2}/shuffleBed -i $*/peaks/final_summit_500.bed \
	-g ${pd}/hg19_canon_chromsizes.txt \
	-excl ${pd}/consensusBlacklist.bed \
	-noOverlapping \
	-seed 12072013 | \
	sed 's/^chr//' > $@; \
	${bedtools}/fastaFromBed -fi ${pd}/hg19.fa -bed $@ -fo $*/peaks/final_summit_500_ctrl.fa

# control peaks final 100bp
make_control_peaks_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/final_summit_100_ctrl.bed)
make_control_peaks : ${make_control_peaks_loop}

%/peaks/final_summit_100_ctrl.bed :
	${bedtools2}/shuffleBed -i $*/peaks/final_summit_100.bed \
	-g ${pd}/hg19_canon_chromsizes.txt \
	-excl ${pd}/consensusBlacklist.bed \
	-noOverlapping \
	-seed 12072013 | \
	sed 's/^chr//' > $@; \
	${bedtools}/fastaFromBed -fi ${pd}/hg19.fa -bed $@ -fo $*/peaks/final_summit_100_ctrl.fa

# control peaks center final 500bp
make_ctrl_center_500_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/final_center_500_ctrl.bed)
make_ctrl_center_500_peaks : ${make_ctrl_center_500_loop}

%/peaks/final_center_500_ctrl.bed :
	${bedtools2}/shuffleBed -i $*/peaks/final_center_500.bed \
	-g ${pd}/hg19_canon_chromsizes.txt \
	-excl ${pd}/consensusBlacklist.bed \
	-noOverlapping \
	-seed 12072013 | \
	sed 's/^chr//' > $@; \
	${bedtools}/fastaFromBed -fi ${pd}/hg19.fa -bed $@ -fo $*/peaks/final_center_500_ctrl.fa

hmod_cov_c_loop = $(foreach lib, ${hmod_libs}, ./merged_peaks/control_padbin_$(lib)_cov.bed)
hmod_cov_c : ${hmod_cov_c_loop}

./merged_peaks/control_padbin_%_cov.bed :
	${bedtools}/coverageBed -abam ${h_chip}/$*/bam/gsc.sorted.bam \
	-b merged_peaks/control_padbin_2kb.bed -counts | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $@

# ------------------------------------------------------------------------------
#  Calculate jaccard index between final CIC peaks and histone mark peaks
# ------------------------------------------------------------------------------
hmod_final_peak_jac_loop = $(foreach lib, ${hmod_libs}, ./jacs/$(lib)_final_peak_jac.tsv)
hmod_final_peak_jac : ${hmod_final_peak_jac_loop}

./jacs/%_final_peak_jac.tsv :
	${bedtools}/bedtools jaccard -a ./merged_peaks/final_peaks_merged.bed \
	-b ${h_chip}/$*/peaks/peaks.sorted.bed | tail -1 > $@

chmm_18_final_peak_jac_loop = $(foreach state, ${chromhmm_18_states}, ./jacs/$(state)/jacs.tsv)
chmm_18_final_peak_jac : ${chmm_18_final_peak_jac_loop}

./jacs/%/jacs.tsv :
	mkdir ./jacs/$*; \
	for i in ${h_chip}/ChromHMM_hg19_18states_output/$*/*; do \
		sample=$$(echo $$i | awk '{gsub("^.*/ChromHMM_hg19_18states_output/$*/", ""); print}' | \
                          awk '{gsub("_18_segments.bed", ""); print}'); \
		jac=$$(${bedtools}/bedtools jaccard -a ./merged_peaks/final_peaks_merged_2.bed -b $$i | tail -1); \
		echo $${sample}	$${jac} >> ./jacs/$*/final_peak_jac.tsv; \
	done

# ------------------------------------------------------------------------------
#  padbin rep2_i_G144 peaks
# ------------------------------------------------------------------------------
A77581/peaks/f_i_G144_padbin2kb.bed :
	cat A77581/peaks/f_i_G144.bed | \
	awk '{print $$1"\t"(($$2 + $$3)/2)}' | \
	awk '{print $$1"\t"($$2 - 2000)"\t"($$2 + 2000)}' | \
	awk '$$2 > 0' | \
	${bedtools}/windowMaker -b - -w 100 -i srcwinnum | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $@

A77581/peaks/f_i_G144_padbin5kb.bed :
	cat A77581/peaks/f_i_G144.bed | \
	awk '{print $$1"\t"(($$2 + $$3)/2)}' | \
	awk '{print $$1"\t"($$2 - 5000)"\t"($$2 + 5000)}' | \
	awk '$$2 > 0' | \
	${bedtools}/windowMaker -b - -w 100 -i srcwinnum | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $@

A77581/peaks/f_i_G144_padbin10kb.bed :
	cat A77581/peaks/f_i_G144.bed | \
	awk '{print $$1"\t"(($$2 + $$3)/2)}' | \
	awk '{print $$1"\t"($$2 - 10000)"\t"($$2 + 10000)}' | \
	awk '$$2 > 0' | \
	${bedtools}/windowMaker -b - -w 100 -i srcwinnum | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $@

A77581/peaks/mq0.01_i_G144_padbin2kb.bed :
	cat A77581/peaks/mq0.01_i_G144.bed | \
	awk '{print $$1"\t"(($$2 + $$3)/2)}' | \
	awk '{print $$1"\t"($$2 - 2000)"\t"($$2 + 2000)}' | \
	awk '$$2 > 0' | \
	${bedtools}/windowMaker -b - -w 100 -i srcwinnum | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $@

A77581/peaks/mq0.01_i_G144_padbin5kb.bed :
	cat A77581/peaks/mq0.01_i_G144.bed | \
	awk '{print $$1"\t"(($$2 + $$3)/2)}' | \
	awk '{print $$1"\t"($$2 - 5000)"\t"($$2 + 5000)}' | \
	awk '$$2 > 0' | \
	${bedtools}/windowMaker -b - -w 100 -i srcwinnum | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $@

A77581/peaks/mq0.01_i_G144_padbin10kb.bed :
	cat A77581/peaks/mq0.01_i_G144.bed | \
	awk '{print $$1"\t"(($$2 + $$3)/2)}' | \
	awk '{print $$1"\t"($$2 - 10000)"\t"($$2 + 10000)}' | \
	awk '$$2 > 0' | \
	${bedtools}/windowMaker -b - -w 100 -i srcwinnum | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $@


# ------------------------------------------------------------------------------
#  Obtain rpkm coverage of IP and input libraries in rep2 peaks overlapping with G144
# ------------------------------------------------------------------------------
./rep2_f_i_G144_rpkm.npz:
	${deeptools}/multiBigwigSummary BED-file -p 10 \
	--bwfiles */bam/gsc.sorted.rpkm.bw \
	--BED A77581/peaks/f_i_G144.bed \
	--labels NHA_rep1_chip NHA_rep1_input NHA_rep2_chip NHA_rep2_input \
	--outRawCounts rep2_f_i_G144_rpkm.tab \
	-out rep2_f_i_G144_rpkm.npz

./rep2_mq0.01_i_G144_rpkm.npz:
	${deeptools}/multiBigwigSummary BED-file -p 10 \
	--bwfiles */bam/gsc.sorted.rpkm.bw \
	--BED A77581/peaks/mq0.01_i_G144.bed \
	--labels NHA_rep1_chip NHA_rep1_input NHA_rep2_chip NHA_rep2_input \
	--outRawCounts rep2_mq0.01_i_G144_rpkm.tab \
	-out rep2_mq0.01_i_G144_rpkm.npz

# ------------------------------------------------------------------------------
#  Obtain rpkm coverage of IP and input libraries across rep2 peaks overlapping with G144
# ------------------------------------------------------------------------------
pads ?= 2kb 5kb 10kb

f_i_G144_padbin_cov_loop = $(foreach pad, ${pads}, A77581/peaks/f_i_G144_$(pad)_rpkm.npz)
f_i_G144_padbin_cov : ${f_i_G144_padbin_cov_loop}

A77581/peaks/f_i_G144_%_rpkm.npz:
	${deeptools}/multiBigwigSummary BED-file -p 20 \
	--bwfiles */bam/gsc.sorted.rpkm.bw \
	--BED A77581/peaks/f_i_G144_padbin$*.bed \
	--labels NHA_rep1_chip NHA_rep1_input NHA_rep2_chip NHA_rep2_input \
	--outRawCounts A77581/peaks/f_i_G144_$*_rpkm.tab \
	-out $@

mq0.01_i_G144_padbin_cov_loop = $(foreach pad, ${pads}, A77581/peaks/mq0.01_i_G144_$(pad)_rpkm.npz)
mq0.01_i_G144_padbin_cov : ${mq0.01_i_G144_padbin_cov_loop}

A77581/peaks/mq0.01_i_G144_%_rpkm.npz:
	${deeptools}/multiBigwigSummary BED-file -p 20 \
	--bwfiles */bam/gsc.sorted.rpkm.bw \
	--BED A77581/peaks/mq0.01_i_G144_padbin$*.bed \
	--labels NHA_rep1_chip NHA_rep1_input NHA_rep2_chip NHA_rep2_input \
	--outRawCounts A77581/peaks/mq0.01_i_G144_$*_rpkm.tab \
	-out $@

# ------------------------------------------------------------------------------
#  Obtain rpkm coverage of hmod ChIP libraries across rep2_i_G144
# ------------------------------------------------------------------------------
f_i_G144_padbin_hmod_cov_loop = $(foreach pad, ${pads}, A77581/peaks/f_i_G144_padbin$(pad)_hmod_rpkm.npz)
f_i_G144_padbin_hmod_cov : ${f_i_G144_padbin_hmod_cov_loop}

A77581/peaks/f_i_G144_padbin%_hmod_rpkm.npz :
	bigwigs=$$(ls ../histone_mod_chip/*/bam/gsc.sorted.rpkm.bw | grep -v input); \
	names=$$(echo $${bigwigs} | awk '{ gsub("/bam/gsc.sorted.rpkm.bw", ""); print }' | \
	awk '{ gsub("../histone_mod_chip/", ""); print }'); \
	${deeptools}/multiBigwigSummary BED-file -p 20 \
	--bwfiles $${bigwigs} \
	--BED A77581/peaks/f_i_G144_padbin$*.bed \
	--labels $${names} \
	--outRawCounts A77581/peaks/f_i_G144_padbin$*_hmod_rpkm.tab \
	-out $@

mq0.01_i_G144_padbin_hmod_cov_loop = $(foreach pad, ${pads}, A77581/peaks/mq0.01_i_G144_padbin$(pad)_hmod_rpkm.npz)
mq0.01_i_G144_padbin_hmod_cov : ${mq0.01_i_G144_padbin_hmod_cov_loop}

A77581/peaks/mq0.01_i_G144_padbin%_hmod_rpkm.npz :
	bigwigs=$$(ls ../histone_mod_chip/*/bam/gsc.sorted.rpkm.bw | grep -v input); \
	names=$$(echo $${bigwigs} | awk '{ gsub("/bam/gsc.sorted.rpkm.bw", ""); print }' | \
	awk '{ gsub("../histone_mod_chip/", ""); print }'); \
	${deeptools}/multiBigwigSummary BED-file -p 20 \
	--bwfiles $${bigwigs} \
	--BED A77581/peaks/mq0.01_i_G144_padbin$*.bed \
	--labels $${names} \
	--outRawCounts A77581/peaks/mq0.01_i_G144_padbin$*_hmod_rpkm.tab \
	-out $@
