# ------------------------------------------------------------------------------
#                           _      _      _
#              		 __(.)< __(.)> __(.)=
#         	         \___)  \___)  \___)
# ------------------------------------------------------------------------------
# the purpose of this pipeline is to perform relevant QC and peak-calling on 
# CIC project histone modification ChIP-seq data.
# 
# it requires as input a file, 'libs.tsv', which has library IDs in the first 
# column and the BAM file location in the second column.
#
# it also requires a second file, 'lib_pairs.tsv', which contains IP and
# input pairs (library IDs) separated by a tab, ie: A77579	A77580
# ------------------------------------------------------------------------------
#              	            _      _      _
#          	    	 __(.)< __(.)> __(.)=
#          	    	 \___)  \___)  \___)
# ------------------------------------------------------------------------------

wd = $(shell pwd)
pd = /projects/marralab_cic_prj/shared_objects/
PERCENT := %
SHELL := /bin/bash

bedtools = /projects/sdlee_prj/installs/bedtools2/bin
bedops = /projects/jtopham_prj/installs/bedops/bin
ChromHMM = /projects/sdlee_prj/installs/ChromHMM
deeptools = /gsc/software/linux-x86_64-centos6/python-2.7.12/bin
finder = /home/mbilenky/bin/Solexa_Java/FindER.1.0.1e.jar
homer = /projects/sdlee_prj/installs/homer/bin
IDR = /gsc/software/linux-x86_64-centos7/idr-2.0.3/bin
java = /gsc/software/linux-x86_64/jre1.8.0_66/bin/java
macs2 = /projects/sdlee_prj/installs/python2.7.12_venv/bin/macs2
picard_dupe = /projects/jtopham_prj/installs/picard-tools-1.114/MarkDuplicates.jar
ssp = /projects/sdlee_prj/installs/SSP/bin/ssp
UCSC = /projects/sdlee_prj/installs/UCSC_tools

# ------------------------------------------------------------------------------
#  set up directory
# ------------------------------------------------------------------------------
libs = $(shell cut -f 1 libs.tsv)
lib_locs = $(shell cat libs.tsv | awk '{print$$1"__"$$2}')
libs_no_input = $(shell cat lib_pairs.tsv | awk '{print $$1}')
samples = $(shell cat samples.tsv)

setupdir = $(foreach lib, ${libs}, $(lib)/bam)
dir : ${setupdir}

# create directory structure
%/bam :
	mkdir $*; mkdir $*/bam;	mkdir $*/peaks;

# import gsc bam files
get_bams = $(foreach lib_loc, ${lib_locs}, $(lib_loc)/bam/gsc.bam)
getbams : ${get_bams}

%/bam/gsc.bam :
	lib=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	loc=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	ln -s $${loc} $${lib}/bam/gsc.bam

# ------------------------------------------------------------------------------
#  process BAM files
# ------------------------------------------------------------------------------
filter_sort = $(foreach lib, ${libs}, $(lib)/bam/gsc.sorted.bam)
sortbams : ${filter_sort}

# filter unmapped reads; sort and index
%/bam/gsc.sorted.bam :
	samtools sort $*/bam/gsc.bam $*/bam/gsc.sorted; \
	samtools index $*/bam/gsc.sorted.bam \
	$*/bam/gsc.sorted.bam.bai

flagstat_loop = $(foreach lib, ${libs}, $(lib)/bam/flagstat.txt)
flagstat : ${flagstat_loop}

# flagstat
%/bam/flagstat.txt :
	samtools flagstat $*/bam/gsc.bam > $@

# filter multi-mapping reads (https://wabi-wiki.scilifelab.se/display/KB/Filter+uniquely+mapped+reads+from+a+BAM+file)
filter_unique_loop = $(foreach lib, ${libs}, $(lib)/bam/gsc.sorted.unique.bam)
filter_unique : ${filter_unique_loop}

%/bam/gsc.sorted.unique.bam :
	samtools view -b -q 5 -F 256 $*/bam/gsc.sorted.bam > $*/bam/tmp.bam; \
	samtools sort $*/bam/tmp.bam $*/bam/gsc.sorted.unique; \
	samtools index $*/bam/gsc.sorted.unique.bam $*/bam/gsc.sorted.unique.bam.bai; \
	rm $*/bam/tmp.bam

rmdupe_loop = $(foreach lib, ${libs}, $(lib)/bam/gsc.nodupe.bam)
rmdupe : ${rmdupe_loop}

%/bam/gsc.nodupe.bam :
	samtools view -bF 4 $*/bam/gsc.sorted.bam > $*/bam/filtered.bam; \
	java -jar ${picard_dupe} I=$*/bam/filtered.bam O=$*/bam/tmp.bam M=$*/bam/dupe.stdout.txt ASSUME_SORTED=true \
	REMOVE_DUPLICATES=true; \
	samtools sort $*/bam/tmp.bam $*/bam/gsc.nodupe; \
	samtools index $@ $@.bai; rm $*/bam/tmp.bam; rm $*/bam/filtered.bam

# ------------------------------------------------------------------------------
#  Check GC Bias
# ------------------------------------------------------------------------------
check_GC = $(foreach lib, ${libs}, $(lib)/bam/GCfreq.txt)
check_GC_bias : ${check_GC}

%/bam/GCfreq.txt :
	computeGCBias -p 20 -b $*/bam/gsc.sorted.bam \
	--effectiveGenomeSize 2747877777 \
	-g /projects/marralab_cic_prj/shared_objects/hg19.2bit \
	-o $@ \
	--biasPlot $*/bam/GCbiasplot.png

# ------------------------------------------------------------------------------
#  Fingerprint
# ------------------------------------------------------------------------------
reps ?= NHA_rep1 NHAA2_rep1 NHA_rep2 NHAA2_rep2 \
	F8_rep1 F8_rep2 F8E10_rep1 F8E10_rep2

fingerprint_loop = $(foreach rep, ${reps}, $(rep)/fingerprints.png)
fingerprint : ${fingerprint_loop}

%/fingerprints.png :
	bams=$$(ls ./$**/bam/gsc.sorted.bam); \
	${deeptools}/plotFingerprint \
	-b $${bams} -p10 \
	--extendReads \
	--labels $$(cat libs.tsv | grep -E $* | awk '{print $$1}' | sed -e "s/^$*_//") \
	--minMappingQuality 30 --skipZeros \
	-T '$* Fingerprints' \
	--plotFile figs/QC/$*_fingerprints.png \
	--outRawCounts objects/$*_fingerprints.tab

# ------------------------------------------------------------------------------
#  Cross Correlation Analysis
# ------------------------------------------------------------------------------
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
#  get coverage of libraries in 175 bp (non-overlapping) bins, in chr 1
# ------------------------------------------------------------------------------
# bin chromosome 1 into 175 bp bins
./chr1_175bp.bed :
	${bedtools}/windowMaker -g ${pd}/hg19_chromsizes.tsv -w 175 | \
	grep ^1 | grep -v ^1[0-9] | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $@

# convert BAM files to BEDs, find read centers, then calculate how many read 
# centers in each 175 bp bin
gcov = $(foreach lib, ${libs}, ./$(lib)/bam/cov_in_genomebins.bed)
gcovs : ${gcov}

%/bam/cov_in_genomebins.bed :
	${bedops}/convert2bed --max-mem=12G --sort-tmpdir=./tmp \
	--input=bam < $*/bam/gsc.bam | cut -f 1-3 | \
	awk '{printf("${PERCENT}i\t${PERCENT}i\n",$$1, ((($$3 - $$2)/2) + \
	$$2))}' | awk '{print $$1"\t"$$2"\t"$$2 + 1}' | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - | \
	${bedtools}/coverageBed -a - -b ./chr1_175bp.bed -counts > $@

# ------------------------------------------------------------------------------
#  merge BAMs across technical replicates
# ------------------------------------------------------------------------------
rep_pairs = $(shell cat rep_pairs.tsv | grep -v "input")

mergebam = $(foreach repp, ${rep_pairs}, ./merged_bam/$(repp).bam)
mergebams : ${mergebam}

./merged_bam/%.bam :
	rep1=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	rep2=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	name=$$(echo $* | awk '{ gsub("_rep.", ""); print }' | \
	awk '{ gsub("__.*", ""); print }'); \
	samtools merge -h $${rep1}/bam/gsc.sorted.bam $${name}/bam/merged.bam \
	$${rep1}/bam/gsc.sorted.bam $${rep2}/bam/gsc.sorted.bam

# ------------------------------------------------------------------------------
#  call peaks using findER
# ------------------------------------------------------------------------------
lib_pairs = $(shell cat lib_pairs.tsv | awk '{print$$1"__"$$2}')

runFinder = $(foreach lib_pair, ${lib_pairs}, $(lib_pair)/peaks.bed)
callpeaks : ${runFinder}

%/peaks.bed :
	ip=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
        input=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	${java} -jar -Xmx120G ${finder} -signalBam \
	./$${ip}/bam/gsc.nodupe.bam -inputBam ./$${input}/bam/gsc.nodupe.bam \
	-out ./$${ip}/peaks/; gunzip -c \
	./$${ip}/peaks/gsc.nodupe.vs.gsc.nodupe.FDR_0.05.FindER.bed.gz | \
        awk '{gsub("chr", "");print}' | grep -v '^GL' | grep -v '^MT' | \
	${bedops}/sort-bed --max-mem 80G --tmpdir ./tmp - \
	> ./$${ip}/peaks/peaks.sorted.bed; \
	cat ./$${ip}/peaks/peaks.sorted.bed | grep -v '^GL' | grep -v '^MT' | \
	./$${ip}/peaks/peaks.sorted.bed; \
	rm ./$${ip}/peaks/gsc.nodupe.vs.gsc.nodupe.FDR_0.05.FindER.bed.gz

# raw coverage in peaks
raw_cov_peaks_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/peaks.sorted.cov.tab)
raw_cov_peaks : ${raw_cov_peaks_loop}

%/peaks/peaks.sorted.cov.tab :
	${deeptools}/multiBamSummary BED-file -p 80 \
	--BED $*/peaks/peaks.sorted.bed \
	--bamfiles $*/bam/gsc.nodupe.bam \
	--outRawCounts $@ \
	--outFileName $*/peaks/peaks.sorted.cov.npz

# RPKM coverage in peaks
rpkm_peaks_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/peaks.rpkm.tab)
rpkm_peaks : ${rpkm_peaks_loop}

%/peaks/peaks.rpkm.tab :
	${deeptools}/multiBigwigSummary BED-file -p 80 \
	--BED $*/peaks/peaks.sorted.rmbl.bed \
	--bwfiles $*/bam/gsc.nodupe.rpkm.bw \
	--outRawCounts $@ \
	--outFileName $*/peaks/peaks.rpkm.npz; \
	awk 'NR>1' $@ | sort -k1,1 -k2,2n > $*/peaks/peaks.rpkm.bed	

input_libs = $(shell cut -f1 libs.tsv | grep "input")

i_rpkm_peaks_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/peaks.i.rpkm.tab)
i_rpkm_peaks : ${i_rpkm_peaks_loop}

%/peaks/peaks.i.rpkm.tab :
	sample=$$(echo $* | awk '{ gsub("h3.*", "input"); print }'); \
	${deeptools}/multiBigwigSummary BED-file -p 80 \
	--BED $*/peaks/peaks.sorted.rmbl.bed \
	--bwfiles $${sample}/bam/gsc.nodupe.rpkm.bw \
	--outRawCounts $@ \
	--outFileName $*/peaks/peaks.i.rpkm.npz; \
	awk 'NR>1' $@ | sort -k1,1 -k2,2n > $*/peaks/peaks.i.rpkm.bed

rpkm_all : rpkm_peaks i_rpkm_peaks

# intersect macs2 peaks across reps
rep_mpeaks_loop = $(foreach pair, ${rep_pairs}, rep_peaks/macs2/$(pair).bed)
rep_mpeaks : ${rep_mpeaks_loop}

rep_peaks/macs2/%.bed :
	rep1=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	rep2=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	name=$$(echo $* | awk '{ gsub("_rep.", ""); print }' | \
	awk '{ gsub("__.*", ""); print }'); \
	${bedtools}/intersectBed -wao -a $${rep1}/peaks/macs2.rmbl.bed \
	-b $${rep2}/peaks/macs2.rmbl.bed > rep_peaks/macs2/$${name}_peaks.rmbl.bed

# overlap and merge finder peaks across reps
m_rep_fold_peaks_loop = $(foreach repp, ${rep_pairs}, merged_peaks/rep_pairs/$(repp).fold.bed)
m_rep_fold_peaks : ${m_rep_fold_peaks_loop}

merged_peaks/rep_pairs/%.fold.bed :
	rep1=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	rep2=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	name=$$(echo $* | awk '{ gsub("_rep.", ""); print }' | \
	awk '{ gsub("__.*", ""); print }'); \
	${bedtools}/intersectBed -wo -a $${rep1}/peaks/peaks.sorted.bed \
	-b $${rep2}/peaks/peaks.sorted.bed > ./merged_peaks/rep_pairs/$${name}_intersect.bed; \
	cat merged_peaks/rep_pairs/$${name}_intersect.bed | awk '{ print$$1"\t"$$2"\t"$$3 }' > tmp/p1; \
	cat merged_peaks/rep_pairs/$${name}_intersect.bed | awk '{ print$$5"\t"$$6"\t"$$7 }' > tmp/p2; \
	cat tmp/p1 tmp/p2 | sort -k1,1 -k2,2n | ${bedtools}/mergeBed -d 1 > merged_peaks/rep_pairs/$${name}.fpeaks.bed

# rmbl, fold filtered (>2) findER peaks
# first overlap peaks between replicates, then merge
finder_fold2_mrep_loop = $(foreach pair, ${rep_pairs}, merged_peaks/rep_pairs/$(pair).fold2.bed)
finder_fold2_mrep : ${finder_fold2_mrep_loop}

merged_peaks/rep_pairs/%.fold2.bed :
	rep1=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	rep2=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	name=$$(echo $* | awk '{ gsub("_rep.", ""); print }' | \
	awk '{ gsub("__.*", ""); print }'); \
	${bedtools}/intersectBed -wo -a $${rep1}/peaks/peaks.fold2.bed \
	-b $${rep2}/peaks/peaks.fold2.bed > ./merged_peaks/rep_pairs/$${name}_intersect.bed; \
	cat merged_peaks/rep_pairs/$${name}_intersect.bed | awk '{ print$$1"\t"$$2"\t"$$3 }' > tmp/p1; \
	cat merged_peaks/rep_pairs/$${name}_intersect.bed | awk '{ print$$5"\t"$$6"\t"$$7 }' > tmp/p2; \
	cat tmp/p1 tmp/p2 | sort -k1,1 -k2,2n | ${bedtools}/mergeBed -d 1 > merged_peaks/rep_pairs/$${name}.fold2.bed

# rmbl, fold filtered (>2) macs2 peaks
# first overlap peaks between replicates, then merge
macs2_fold2_mrep_loop = $(foreach pair, ${rep_pairs}, merged_peaks/rep_pairs/macs2/$(pair).fold2.bed)
macs2_fold2_mrep : ${macs2_fold2_mrep_loop}

merged_peaks/rep_pairs/macs2/%.fold2.bed :
	rep1=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	rep2=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	name=$$(echo $* | awk '{ gsub("_rep.", ""); print }' | \
	awk '{ gsub("__.*", ""); print }'); \
	${bedtools}/intersectBed -wo -a $${rep1}/peaks/macs2.rmbl.fold2.bed \
	-b $${rep2}/peaks/macs2.rmbl.fold2.bed > ./merged_peaks/rep_pairs/macs2/$${name}_intersect.bed; \
	cat merged_peaks/rep_pairs/macs2/$${name}_intersect.bed | awk '{ print$$1"\t"$$2"\t"$$3 }' > tmp/p1; \
	cat merged_peaks/rep_pairs/macs2/$${name}_intersect.bed	| awk '{ print$$5"\t"$$6"\t"$$7 }' > tmp/p2; \
	cat tmp/p1 tmp/p2 | sort -k1,1 -k2,2n | ${bedtools}/mergeBed -d 1 > merged_peaks/rep_pairs/macs2/$${name}.rmbl.fold2.bed

# Annotate peaks
ips = $(shell cat lib_pairs.tsv | awk '{print$$1}')

anno_peak = $(foreach ip, ${ips}, $(ip)/peaks/peaks.sorted.anno.bed)
anno_peaks : ${anno_peak}

%/peaks/peaks.sorted.anno.bed :
	${bedtools}/closestBed -a $*/peaks/peaks.sorted.bed -b ${pd}/Homo_sapiens_GRCh37_genes.sorted.bed | \
	awk '{print $$1"\t"$$2"\t"$$3"\t"$$11"\t"$$4"\t"$$10}' > $@

# ------------------------------------------------------------------------------
#  filter blacklist
# ------------------------------------------------------------------------------
m_samples = $(shell cat all_samples_repmerged.tsv)

rm_blacklist_loop = $(foreach sp, ${m_samples}, rep_peaks/$(sp).rmbl.bed)
rm_blacklist : ${rm_blacklist_loop}

# cat ../shared_objects/consensusBlacklist.bed | sed 's/^chr//' > tmp/bl.bed
rep_peaks/%.rmbl.bed :
	${bedtools}/subtractBed -A \
	-a rep_peaks/$*_peaks.bed \
	-b tmp/bl.bed > rep_peaks/$*_peaks.rmbl.bed

# remove bl from macs2 peaks
rmbl_macs2_peaks_loop = $(foreach ip, ${ips}, $(ip)/peaks/macs2_peaks.rmbl.bed)
rmbl_macs2_peaks : ${rmbl_macs2_peaks_loop}

%/peaks/macs2_peaks.rmbl.bed :
	${bedtools}/subtractBed -A \
	-a $*/peaks/macs2_q0.01_peaks.narrowPeak \
	-b tmp/bl.bed > $*/peaks/tmp.bed; \
	sort -k1,1 -k2,2n $*/peaks/tmp.bed > $@; rm $*/peaks/tmp.bed

# remove bl from peaks.sorted.tab
rm_bl_fpeaks_loop = $(foreach ip, ${ips}, $(ip)/peaks/peaks.sorted.rmbl.bed)
rm_bl_fpeaks : ${rm_bl_fpeaks_loop}

%/peaks/peaks.sorted.rmbl.bed :
	${bedtools}/subtractBed -A \
	-a $*/peaks/peaks.sorted.bed \
	-b tmp/bl.bed > $*/peaks/tmp.bed; \
	sort -k1,1 -k2,2n $*/peaks/tmp.bed > $@; rm $*/peaks/tmp.bed

# remove bl from peaks.cov.tab
rm_bl_loop = $(foreach ip, ${ips}, $(ip)/peaks/peaks.sorted.cov.rmbl.bed)
rm_bl : ${rm_bl_loop}

%/peaks/peaks.sorted.cov.rmbl.bed :
	sort -k1,1 -k2,2n $*/peaks/peaks.sorted.cov.tab > $*/peaks/peaks.sorted.cov.bed; \
	${bedtools}/subtractBed -A \
	-a $*/peaks/peaks.sorted.cov.bed \
	-b tmp/bl.bed > $*/peaks/tmp.bed; \
	sort -k1,1 -k2,2n $*/peaks/tmp.bed > $@; rm $*/peaks/tmp.bed

# ------------------------------------------------------------------------------
#  IDR analysis
# ------------------------------------------------------------------------------
broad_rep_pairs ?= $(shell cat rep_pairs.tsv | awk '{print$$1__$$2}' | grep -v "h3k27ac" | grep -v "h3k4" | grep -v "input")
narrow_rep_pairs ?= $(shell cat rep_pairs.tsv | awk '{print$$1__$$2}' | grep -v "h3k27me3" | grep -v "h3k36" | grep -v "h3k9" | grep -v "input")

idr_broad_loop = $(foreach pair, ${broad_rep_pairs}, IDR/$(pair).broad.idr.bed)
idr_broad : ${idr_broad_loop}

IDR/%.broad.idr.bed :
	rep1=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	rep2=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	name=$$(echo $* | awk '{ gsub("_rep.", ""); print }' | \
	awk '{ gsub("__.*", ""); print }'); \
	idr --sample $${rep1}/peaks/macs2.bed $${rep2}/peaks/macs2.bed \
	--input-file-type broadPeak \
	--output-file-type broadPeak \
	--output-file IDR/$${name}.idr.bed \
	--max-iter 300 \
	--log-output-file IDR/$${name}.idr.log

idr_narrow_loop = $(foreach pair, ${narrow_rep_pairs}, IDR/$(pair).narrow.idr.bed)
idr_narrow : ${idr_narrow_loop}

IDR/%.narrow.idr.bed :
	rep1=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	rep2=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	name=$$(echo $* | awk '{ gsub("_rep.", ""); print }' | \
	awk '{ gsub("__.*", ""); print }'); \
	idr --sample $${rep1}/peaks/macs2.bed $${rep2}/peaks/macs2.bed \
	--input-file-type narrowPeak \
	--output-file-type narrowPeak \
	--output-file IDR/$${name}.idr.bed \
	--max-iter 300 \
	--log-output-file IDR/$${name}.idr.log

macs2_idr : idr_broad idr_narrow

# ------------------------------------------------------------------------------
#  raw coverage in peaks (rmbl, fold>1)
# ------------------------------------------------------------------------------
cov_peaks_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/peaks.fold.filt.cov.tab)
cov_peaks : ${cov_peaks_loop}

%/peaks/peaks.fold.filt.cov.tab :
	${deeptools}/multiBamSummary BED-file -p 120 \
	--BED $*/peaks/peaks.fold.filt.bed \
	--bamfiles $*/bam/gsc.nodupe.bam \
	--outRawCounts $@ \
	--outFileName $*/peaks/peaks.fold.filt.cov.npz; \
	awk 'NR>1' $@ | sort -k1,1 -k2,2n > $*/peaks/peaks.fold.filt.cov.bed

cov_fold2_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/macs2.rmbl.fold2.cov.tab)
cov_fold2 : ${cov_fold2_loop}

%/peaks/macs2.rmbl.fold2.cov.tab :
	${deeptools}/multiBamSummary BED-file -p 80 \
	--BED $*/peaks/macs2.rmbl.fold2.bed \
	--bamfiles $*/bam/gsc.nodupe.bam \
	--outRawCounts $@ \
	--outFileName $*/peaks/macs2.rmbl.fold2.cov.npz; \
	awk 'NR>1' $@ | sort -k1,1 -k2,2n > $*/peaks/macs2.rmbl.fold2.cov.bed

# ------------------------------------------------------------------------------
#  Union of all FindER peaks for each mark
# ------------------------------------------------------------------------------
hmods ?= h3k4me1 h3k4me3 h3k27ac h3k27me3 h3k36me3 h3k9me3

union_finder_loop = $(foreach hmod, ${hmods}, finder_mpeaks/$(hmod)/union.bed)
union_finder : ${union_finder_loop}

./finder_mpeaks/%/union.bed :
	mkdir ./finder_mpeaks/$* ; \
	cat ./*$*/peaks/peaks.sorted.bed | sort -k1,1 -k2,2n | \
	${bedtools}/mergeBed > $@

u_fpeaks_loop = $(foreach hmod, ${hmods}, merged_peaks/$(hmod)/union.rmbl.bed)
u_fpeaks : ${u_fpeaks_loop}

merged_peaks/%/union.rmbl.bed :
	cat ./merged_peaks/rep_pairs/*$*.rmbl.bed | sort -k1,1 -k2,2n | \
	${bedtools}/mergeBed > $@

# union of all peaks in all samples
u_fpeaks_fold_loop = $(foreach hmod, ${hmods}, merged_peaks/$(hmod)/u_fpeaks_f2.bed)
u_fpeaks_fold : ${u_fpeaks_fold_loop}

merged_peaks/%/u_fpeaks_f2.bed :
	cat *$*/peaks/peaks.fold2.bed | sort -k1,1 -k2,2n | \
	${bedtools}/mergeBed > $@

# union	of rep merged peaks
u_fpeaks_fold_loop = $(foreach hmod, ${hmods}, merged_peaks/$(hmod)/union.finder.fold2.bed)
u_fpeaks_fold : ${u_fpeaks_fold_loop}

merged_peaks/%/union.finder.fold2.bed :
	cat ./merged_peaks/rep_pairs/*$*.fold2.bed | sort -k1,1 -k2,2n | \
	${bedtools}/mergeBed > $@

# ------------------------------------------------------------------------------
#  SPKM in union FindER peaks (rmbl, FE > 2), keep replicates separate
# ------------------------------------------------------------------------------
u_fpeaks_f2_spkm_loop = $(foreach hmod, ${hmods}, merged_peaks/$(hmod)/u_fpeaks_f2.spkm.tab)
u_fpeaks_f2_spkm : ${u_fpeaks_f2_spkm_loop}

merged_peaks/%/u_fpeaks_f2.spkm.tab :
	bigwigs=$$(ls *$*/bam/gsc.nodupe.rpkm.bw); \
	names=$$(echo $${bigwigs} | awk '{ gsub("$*/bam/gsc.nodupe.rpkm.bw", ""); print}'); \
	${deeptools}/multiBigwigSummary BED-file -p 60 \
	--bwfiles $${bigwigs} \
	--BED merged_peaks/$*/u_fpeaks_f2.bed \
	--labels $${names} \
	--outRawCounts $@ \
	-out merged_peaks/$*/u_fpeaks_f2.spkm.npz


# union of rep merged macs2 peaks (rmbl, fold > 2)
u_mpeaks_fold2_loop = $(foreach hmod, ${hmods}, merged_peaks/$(hmod)/union.macs2.rmbl.fold2.bed)
u_mpeaks_fold2 : ${u_mpeaks_fold2_loop}

merged_peaks/%/union.macs2.rmbl.fold2.bed :
	cat ./merged_peaks/rep_pairs/macs2/*$*.rmbl.fold2.bed | sort -k1,1 -k2,2n | \
	${bedtools}/mergeBed > $@

# ------------------------------------------------------------------------------
#  presence of union peak in each sample
# ------------------------------------------------------------------------------
peak_spec_loop = $(foreach mod, ${hmods}, ./merged_peaks/$(mod)/union.finder.fold2.spec.bed)
peak_spec : ${peak_spec_loop}

./merged_peaks/%/union.finder.fold2.spec.bed :
	${bedtools}/intersectBed -loj -filenames \
	-a ./merged_peaks/$*/union.finder.fold2.bed \
	-b ./merged_peaks/rep_pairs/*$*.fold2.bed > $@

# ------------------------------------------------------------------------------
#  Coverage within union FindER Peaks
# ------------------------------------------------------------------------------
u_finder_cov_loop = $(foreach hmod, ${hmods}, finder_mpeaks/$(hmod)/union_cov.tab)
u_finder_cov : ${u_finder_cov_loop}

./finder_mpeaks/%/union_cov.tab :
	bams=$$(ls *$*/bam/gsc.sorted.bam); \
	names=$$(echo $${bams} | awk '{ gsub("_$*/bam/gsc.sorted.bam", ""); print }'); \
	${deeptools}/multiBamSummary BED-file -p 30 \
	--bamfiles $${bams} \
	--BED finder_mpeaks/$*/union.bed \
	--labels $${names} \
	--outRawCounts $@ \
	-out finder_mpeaks/$*/union_cov.npz

u_fpeaks_cov_loop = $(foreach hmod, ${hmods}, merged_peaks/$(hmod)/union.rmbl.tab)
u_fpeaks_cov : ${u_fpeaks_cov_loop}

./merged_peaks/%/union.rmbl.tab :
	bams=$$(ls *$*/bam/gsc.nodupe.bam); \
	names=$$(echo $${bams} | awk '{ gsub("_$*/bam/gsc.nodupe.bam", ""); print }'); \
	${deeptools}/multiBamSummary BED-file -p 80 \
	--bamfiles $${bams} \
	--BED merged_peaks/$*/union.rmbl.bed \
	--labels $${names} \
	--outRawCounts $@ \
	-out merged_peaks/$*/union.rmbl.npz

# rpkm within union FindER peaks (using merged bigwigs between reps)
hmods_input = $(shell echo ${hmods} input)

u_f_rpkm_loop = $(foreach mod, ${hmods_input}, merged_peaks/$(mod)/union.finder.fold2.rpkm.tab)
u_f_rpkm : ${u_f_rpkm_loop}

merged_peaks/%/union.finder.fold2.rpkm.tab :
	bigwigs=$$(ls merged_bw/rpkm/$*/*.rpkm.bw); \
	names=$$(echo $${bigwigs} | awk '{ gsub("merged_bw/rpkm/$*/", ""); print}' | \
	awk '{ gsub(".rpkm.bw", ""); print }'); \
	${deeptools}/multiBigwigSummary BED-file -p 80 \
	--bwfiles $${bigwigs} \
	--BED merged_peaks/$*/union.finder.fold2.bed \
	--labels $${names} \
	--outRawCounts $@ \
	-out merged_peaks/$*/union.finder.fold2.rpkm.npz; \
	awk 'NR>1' | sort -k1,1 -k2,2n > $@

# input rpkm within union FindER peaks
input_u_f_rpkm_loop = $(foreach hmod, ${hmods}, merged_peaks/$(hmod)/union.finder.fold2.i.rpkm.tab)
input_u_f_rpkm : ${input_u_f_rpkm_loop}

merged_peaks/%/union.finder.fold2.i.rpkm.tab :
	bigwigs=$$(ls merged_bw/rpkm/input/*.rpkm.bw); \
	names=$$(echo $${bigwigs} | awk '{ gsub("merged_bw/rpkm/input/", ""); print}' | \
	awk '{ gsub(".rpkm.bw", ""); print }'); \
	${deeptools}/multiBigwigSummary BED-file -p 80 \
	--bwfiles $${bigwigs} \
	--BED merged_peaks/$*/union.finder.fold2.bed \
	--labels $${names} \
	--outRawCounts $@ \
	-out merged_peaks/$*/union.finder.fold2.i.rpkm.npz; \
	awk 'NR>1' $@ | sort -k1,1 -k2,2n > merged_peaks/$*/union.finder.fold2.i.rpkm.bed

# coverage within peaks (rmbl, fold>2, merged_rep) for DE analysis
u_fpeaks_fold2_cov_loop = $(foreach hmod, ${hmods}, merged_peaks/$(hmod)/union.finder.fold2.tab)
u_fpeaks_fold2_cov : ${u_fpeaks_fold2_cov_loop}

./merged_peaks/%/union.finder.fold2.tab :
	bams=$$(ls *$*/bam/gsc.nodupe.bam); \
	names=$$(echo $${bams} | awk '{ gsub("_$*/bam/gsc.nodupe.bam", ""); print }'); \
	${deeptools}/multiBamSummary BED-file -p 60 \
	--bamfiles $${bams} \
	--BED merged_peaks/$*/union.finder.fold2.bed \
	--labels $${names} \
	--outRawCounts $@ \
	-out merged_peaks/$*/union.finder.fold2.npz

sort_u_f_rpkm_loop = $(foreach hmod, ${hmods}, merged_peaks/$(hmod)/union.rmbl.i.rpkm.bed)
sort_u_f_rpkm : ${sort_u_f_rpkm_loop}

merged_peaks/%/union.rmbl.i.rpkm.bed :
	awk 'NR>1' merged_peaks/$*/union.rmbl.i.rpkm.tab | sort -k1,1 -k2,2n > $@

# log2 coverage (relative to input) for union FindER peaks (using merged bigwigs between reps)
u_f_log2_loop = $(foreach hmod, ${hmods}, merged_peaks/$(hmod)/union.rmbl.log2.tab)
u_f_log2 : ${u_f_log2_loop}

merged_peaks/%/union.rmbl.log2.tab :
	bigwigs=$$(ls merged_bw/input_norm/$*/*nodupe.log2.bw); \
	names=$$(echo $${bigwigs} | awk '{ gsub("merged_bw/input_norm/$*/", ""); print}' | \
	awk '{ gsub(".nodupe.log2.bw", ""); print }'); \
	${deeptools}/multiBigwigSummary BED-file -p 100 \
	--bwfiles $${bigwigs} \
	--BED merged_peaks/$*/union.rmbl.bed \
	--labels $${names} \
	--outRawCounts $@ \
	-out merged_peaks/$*/union.rmbl.log2.npz

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

broad_lib_pairs ?= $(shell cat lib_pairs.tsv | awk '{print$$1"__"$$2}' | grep -v "h3k27ac" | grep -v "h3k4")

macs2_broad_loop = $(foreach pair, ${broad_lib_pairs}, $(pair)/macs2.broadPeak)
macs2_broad : ${macs2_broad_loop}

%/macs2.broadPeak :
	ip=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	input=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	${macs2} callpeak -B \
	-t $${ip}/bam/gsc.sorted.bam \
	-c $${input}/bam/gsc.sorted.bam \
	--broad \
	--broad-cutoff 0.05 \
	-f BAMPE \
	-g 2736124973 \
	-n broad_mq0.05 \
	--outdir $${ip}/peaks

narrow_libs ?= $(shell cut -f1 libs.tsv | grep -v "h3k27me3" | grep -v "h3k9" | grep -v "h3k36me3" | grep -v "input")
broad_libs ?= $(shell cut -f1 libs.tsv | grep -v "h3k27ac" | grep -v "h3k4" | grep -v "input")

# cp narrow and broad macs2 peak files for consistent naming
macs2_broad_cp_loop = $(foreach lib, ${broad_libs}, $(lib)/peaks/broad_macs2.bed)
macs2_broad_cp : ${macs2_broad_cp_loop}

%/peaks/broad_macs2.bed :
	cp $*/peaks/broad_mq0.05.broadPeak $*/peaks/macs2.bed

macs2_narrow_cp_loop = $(foreach lib, ${narrow_libs}, $(lib)/peaks/narrow_macs2.bed)
macs2_narrow_cp : ${macs2_narrow_cp_loop}

%/peaks/narrow_macs2.bed :
	cp $*/peaks/macs2_q0.01_peaks.narrowPeak $*/peaks/macs2.bed

# remove bl from macs2 peaks
rmbl_macs2_loop = $(foreach ip, ${ips}, $(ip)/peaks/macs2.rmbl.bed)
rmbl_macs2 : ${rmbl_macs2_loop}

%/peaks/macs2.rmbl.bed :
	${bedtools}/subtractBed -A \
	-a $*/peaks/macs2.bed \
	-b tmp/bl.bed > $*/peaks/tmp.bed; \
	sort -k1,1 -k2,2n $*/peaks/tmp.bed > $@; rm $*/peaks/tmp.bed

macs2_cp : macs2_broad_cp macs2_narrow_cp rmbl_macs2

# Raw coverage in macs2.rmbl.bed
raw_cov_macs2_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/macs2.rmbl.cov.tab)
raw_cov_macs2 : ${raw_cov_macs2_loop}

%/peaks/macs2.rmbl.cov.tab :
	${deeptools}/multiBamSummary BED-file -p 100 \
	--BED $*/peaks/macs2.rmbl.bed \
	--bamfiles $*/bam/gsc.nodupe.bam \
	--outRawCounts $@ \
	--outFileName $*/peaks/macs2.rmbl.cov.npz; \
	awk 'NR>1' $@ | sort -k1,1 -k2,2n > $*/peaks/macs2.rmbl.cov.bed

# RPKM in macs2.rmbl.bed
rpkm_macs2_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/macs2.rmbl.rpkm.tab)
rpkm_macs2 : ${rpkm_macs2_loop}

%/peaks/macs2.rmbl.rpkm.tab :
	${deeptools}/multiBigwigSummary BED-file -p 60 \
	--BED $*/peaks/macs2.rmbl.bed \
	--bwfiles $*/bam/gsc.nodupe.rpkm.bw \
	--outRawCounts $@ \
	--outFileName $*/peaks/macs2.rmbl.rpkm.npz; \
	awk 'NR>1' $@ | sort -k1,1 -k2,2n > $*/peaks/macs2.rmbl.rpkm.bed

# input RPKM in macs2.rmbl.bed
i_rpkm_macs2_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/macs2.rmbl.i.rpkm.tab)
i_rpkm_macs2 : ${i_rpkm_macs2_loop}

%/peaks/macs2.rmbl.i.rpkm.tab :
	sample=$$(echo $* | awk '{ gsub("h3.*", "input"); print }'); \
	${deeptools}/multiBigwigSummary BED-file -p 60 \
	--BED $*/peaks/macs2.rmbl.bed \
	--bwfiles $${sample}/bam/gsc.nodupe.rpkm.bw \
	--outRawCounts $@ \
	--outFileName $*/peaks/macs2.rmbl.i.rpkm.npz; \
	awk 'NR>1' $@ | sort -k1,1 -k2,2n > $*/peaks/macs2.rmbl.i.rpkm.bed

rpkm_macs2_rmbl : rpkm_macs2 i_rpkm_macs2

# ------------------------------------------------------------------------------
#  Union of all MACS2 peaks for each mark
# ------------------------------------------------------------------------------
narrow_libs ?= $(shell cut -f1 libs.tsv | grep -v "h3k27me3" | grep -v "h3k9" | grep -v "h3k36me3" | grep -v "input")
broad_libs ?= $(shell cut -f1 libs.tsv | grep -v "h3k27ac" | grep -v "h3k4" | grep -v "input")

narrow_marks ?= h3k4me1 h3k4me3 h3k27ac
broad_marks ?= h3k27me3 h3k36me3 h3k9me3

union_narrow_loop = $(foreach mark, ${narrow_marks}, macs2_mpeaks/narrow/$(mark))
union_narrow : ${union_narrow_loop}

./macs2_mpeaks/narrow/% :
	mkdir ./macs2_mpeaks/$* ; \
	cat ./*$*/peaks/macs2_q0.01_peaks.narrowPeak | sort -k1,1 -k2,2n | \
	${bedtools}/mergeBed > ./macs2_mpeaks/$*/union.bed
	
union_broad_loop = $(foreach mark, ${broad_marks}, macs2_mpeaks/broad/$(mark))
union_broad : ${union_broad_loop}

./macs2_mpeaks/broad/% :
	mkdir ./macs2_mpeaks/$* ; \
	cat ./*$*/peaks/broad_mq0.05.broadPeak | sort -k1,1 -k2,2n | \
	${bedtools}/mergeBed > ./macs2_mpeaks/$*/union.bed

macs2_union : union_narrow union_broad

# ------------------------------------------------------------------------------
#  Coverage within union MACS2 Peaks
# ------------------------------------------------------------------------------
u_macs2_cov_loop = $(foreach hmod, ${hmods}, macs2_mpeaks/$(hmod)/union_cov.tab)
u_macs2_cov : ${u_macs2_cov_loop}

macs2_mpeaks/%/union_cov.tab :
	bams=$$(ls *$*/bam/gsc.sorted.bam); \
	names=$$(echo $${bams} | awk '{ gsub("_$*/bam/gsc.sorted.bam", ""); print }'); \
	${deeptools}/multiBamSummary BED-file -p 40 \
	--bamfiles $${bams} \
	--BED macs2_mpeaks/$*/union.bed \
	--labels $${names} \
	--outRawCounts macs2_mpeaks/$*/union_cov.tab \
	-out macs2_mpeaks/$*/union_cov.npz

# ------------------------------------------------------------------------------
#  coverage within peaks
# ------------------------------------------------------------------------------
peakcov = $(foreach lib_pair, ${lib_pairs}, $(lib_pair)/peaks/cov_in_peaks.bed)
peakcovs : ${peakcov}

# raw coverage in peaks
%/peaks/cov_in_peaks.bed :
	ip=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	input=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	${bedtools}/coverageBed -abam $${ip}/bam/gsc.bam -b \
	$${ip}/peaks/peaks.sorted.bed -counts | \
	${bedops}/sort-bed --max-mem 10G --tmpdir ./tmp - > \
	$${ip}/bam/cov_in_$${ip}_peaks.bed; \
	${bedtools}/coverageBed -abam $${input}/bam/gsc.bam -b \
        $${ip}/peaks/peaks.sorted.bed -counts | \
        ${bedops}/sort-bed --max-mem 10G --tmpdir ./tmp - > \
	$${input}/bam/cov_in_$${ip}_peaks.bed

peakcov_mq0.01_loop = $(foreach lib_pair, ${lib_pairs}, $(lib_pair)/cov_in_mq0.01_peaks.bed)
peakcov_mq0.01 : ${peakcov_mq0.01_loop}

%/cov_in_mq0.01_peaks.bed :
	ip=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	input=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	${bedtools}/coverageBed -abam $${ip}/bam/gsc.sorted.bam -b \
	$${ip}/peaks/macs2_q0.01_peaks.narrowPeak -counts | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $${ip}/peaks/cov_in_mq0.01_peaks.bed; \
	${bedtools}/coverageBed -abam $${input}/bam/gsc.sorted.bam -b \
	$${ip}/peaks/macs2_q0.01_peaks.narrowPeak -counts | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $${input}/peaks/cov_in_mq0.01_peaks.bed; \
	paste $${ip}/peaks/cov_in_mq0.01_peaks.bed <(awk '{print $$5}' $${input}/peaks/cov_in_mq0.01_peaks.bed) > \
	$${ip}/peaks/both_cov_in_mq0.01_peaks.bed; \
	rm $${input}/peaks/cov_in_mq0.01_peaks.bed; rm $${ip}/peaks/cov_in_mq0.01_peaks.bed 

padbin_mq0.01_peaks_loop = $(foreach lib, ${libs_no_input}, $(lib)/peaks/macs2_q0.01_peaks_padbin.bed)
padbin_mq0.01_peaks : ${padbin_mq0.01_peaks_loop}

%/peaks/macs2_q0.01_peaks_padbin.bed :
	cat $*/peaks/macs2_q0.01_peaks.narrowPeak | \
	awk '{print $$1"\t"int(($$2 + $$3)/2)}' | \
	awk '{print $$1"\t"($$2 - 2000)"\t"($$2 + 2000)}' | \
	awk '$$2 > 0' | \
	${bedtools}/windowMaker -b - -w 100 -i srcwinnum | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $@

cov_around_macs2_peaks_loop = $(foreach lib_pair, ${lib_pairs}, $(lib_pair)/cov_around_macs2_peaks.bed)
cov_around_macs2_peaks : ${cov_around_macs2_peaks_loop}

%/cov_around_macs2_peaks.bed :
	ip=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	input=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	${bedtools}/coverageBed -abam $${ip}/bam/gsc.sorted.bam \
	-b $${ip}/peaks/macs2_q0.01_peaks_padbin.bed -counts | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $${ip}/peaks/macs2_q0.01_peaks_padbin_cov.bed; \
	${bedtools}/coverageBed -abam $${input}/bam/gsc.sorted.bam \
	-b $${ip}/peaks/macs2_q0.01_peaks_padbin.bed -counts | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $${input}/peaks/macs2_q0.01_peaks_padbin_cov.bed; \
	paste $${ip}/peaks/macs2_q0.01_peaks_padbin_cov.bed <(awk '{print $$5}' $${input}/peaks/macs2_q0.01_peaks_padbin_cov.bed) > \
	$${ip}/peaks/both_macs2_q0.01_peaks_padbin_cov.bed

# ------------------------------------------------------------------------------
#  TSS meta plots
# ------------------------------------------------------------------------------
# pad and bin TSS regions
./tss_padbin.bed :
	cat ${pd}/tss_pad_1500_500.bed | \
	awk '{printf("${PERCENT}s\t${PERCENT}i\t${PERCENT}i\t${PERCENT}s\n", \
	$$1, ($$10 - 5000), ($$10 + 5000), $$6)}' | awk '$$2 > 0' | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - | \
	${bedtools}/windowMaker -b - -w 100 -i srcwinnum | \
        ${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $@

# coverage of IPs in TSS bins
tsscov = $(foreach lib_pair, ${lib_pairs}, $(lib_pair)/peaks/cov_across_tss.bed)
tsscovs : ${tsscov}

%/peaks/cov_across_tss.bed :
	ip=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
        ${bedtools}/coverageBed -abam ./$${ip}/bam/gsc.sorted.bam \
	-b ./tss_padbin.bed -counts | \
	tr -s '_' $$'\t' > ./$${ip}/peaks/cov_across_tss.bed

# coverage of IPs in TSS bins using deeptools
cov_tss_padbin_5kb.npz :
	bams=$$(ls */bam/gsc.sorted.bam); \
	names=$$(echo $${bams} | awk '{ gsub("/bam/gsc.sorted.bam", ""); print }'); \
	${deeptools}/multiBamSummary BED-file -p 40 \
	-b $${bams} \
	--BED tss_padbin.bed \
	--labels $${names} \
	--outRawCounts cov_tss_padbin_5kb.tab \
	-out $@

# coverage of IPs in TSS bins no dupes
# bams=$$(ls */bam/gsc.nodupe.bam); \
# names=$$(echo $${bams} | awk '{ gsub("/bam/gsc.nodupe.bam", ""); print }'); \

cov_tss_padbin_5kb_nodup.npz :
	${deeptools}/multiBamSummary BED-file -p 120 \
	-b NHA_rep1_h3k27ac/bam/gsc.nodupe.bam NHA_rep2_h3k27ac/bam/gsc.nodupe.bam \
	--BED tss_padbin_canon.bed \
	--labels rep1 rep2 \
	--outRawCounts cov_tss_padbin_5kb_nodup.tab \
	--transcript_id_designator \
	-out $@


cov_tss_padbin_5kb_rpkm.npz :
	bws=$$(ls */bam/gsc.sorted.rpkm.bw); \
	names=$$(echo $${bws} | awk '{ gsub("/bam/gsc.sorted.rpkm.bw", ""); print }'); \
	${deeptools}/multiBigwigSummary BED-file -p 40 \
	-b $${bws} \
	--BED tss_padbin_canon.bed \
	--labels $${names} \
	--outRawCounts cov_tss_padbin_5kb_rpkm.tab \
	-out $@

cov_tss_padbin_5kb_meanrep_rpkm.npz :
	bws=$$(ls merged_bw/*.bw); \
	names=$$(echo $${bws} | awk '{ gsub("merged_bw/", ""); print }' | \
	awk '{ gsub(".rpkm.bw", ""); print }'); \
	echo $${names}; \
	${deeptools}/multiBigwigSummary BED-file -p 80 \
	-b $${bws} \
	--BED tss_padbin_canon.bed \
	--labels $${names} \
	--outRawCounts cov_tss_padbin_5kb_rpkm.tab \
	-out $@

# coverage of IPs in CIC positive control genes
poscov = $(foreach ip, ${ips}, $(ip)/peaks/cov_across_pos_genes.bed)
poscovs : ${poscov}

%/peaks/cov_across_pos_genes.bed :
	awk '{if ($$4 == "ETV1" || $$4  == "ETV4" || $$4  == "ETV5" || $$4  == "DUSP4" || $$4  == "GPR3" \
	|| $$4  == "SPRY4" || $$4  == "SHC3" || $$4  == "SHC4") {print}}' $*/peaks/cov_across_tss.bed > $@

# ------------------------------------------------------------------------------
#  jaccard indeces between lists of peaks for each IP
# ------------------------------------------------------------------------------
jac_pairs = $(shell libs=$$(cut -f 1 ./libs.tsv); for i in $$libs; \
  	      do for j in $$libs; do echo $$i"__"$$j; done; done | grep -v input)

jac = $(foreach jp, ${jac_pairs}, ./jaccards/$(jp).tsv)
jacs : ${jac}

./jaccards/%.tsv :
	p1=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	p2=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	${bedtools}/bedtools jaccard -a ./$${p1}/peaks/peaks.sorted.bed \
	-b ./$${p2}/peaks/peaks.sorted.bed | tail -1 > $@

# ------------------------------------------------------------------------------
#  Bam to bedGraph
# ------------------------------------------------------------------------------
bedGraph = $(foreach lib, ${libs}, $(lib)/bam/bam.bedGraph)
sortbedGraph = $(foreach lib, ${libs}, $(lib)/bam/bam.sorted.bedGraph)
BigWig = $(foreach lib, ${libs}, $(lib)/bam/bam.sorted.bw)

bedGraphs : ${bedGraph}
sortbedGraphs : ${sortbedGraph}
BigWigs : ${BigWig}

%/bam/bam.bedGraph :
	${bedtools}/genomeCoverageBed -ibam $*/bam/gsc.sorted.bam -bg  > $@

%/bam/bam.sorted.bedGraph :
	grep -v '^MT' $*/bam/bam.bedGraph | grep -v '^GL' | awk '{print "chr"$$0}' | sort -k1,1 -k2,2n > $@

%/bam/bam.sorted.bw :
	/projects/echun_prj/RT/bin/UCSC_Tools/bedGraphToBigWig $*/bam/bam.sorted.bedGraph /projects/echun_prj/RT/bin/UCSC_Tools/hg19.chrom.sizes $@; \
	cp $@ /gsc/www/bcgsc.ca/downloads/sdlee/CIC/$*.bam.sorted.bw; \
	echo "track type=bigWig name='$*' description='$*' \
	color=0,0,255 bigDataUrl=http://bigwigviewer:bwvbwv@www.bcgsc.ca/downloads/sdlee/CIC/$*.bam.sorted.bw" >> upload

cleanbg :
	for dir in $(libs) ; do \
		rm $$dir/bam/bam.bedGraph; \
	done

cleanbw :
	for dir in $(libs) ; do \
		rm $$dir/bam/bam.sorted.bw; \
	done

# ------------------------------------------------------------------------------
#  ChromHMM
# ------------------------------------------------------------------------------
ChromHMMdir :
	mkdir ChromHMMdir2; \
	
	for dir in ${libs}; do \
		${bedtools}/bamToBed -i $$dir/bam/gsc.nodupe.bam | \
		awk '{print "chr"$$0}' > ChromHMMdir2/$$dir.bed; \
	done

binarizeBeds :
	mkdir binarized_beds2; \

	${java} -jar -Xmx100G ${ChromHMM}/ChromHMM.jar BinarizeBed /projects/sdlee_prj/installs/ChromHMM/CHROMSIZES/hg19.txt \
	ChromHMMdir2 cellmarkfiletable.txt binarized_beds2

# run ChromHMM LearnModel to call 15 chromatin states
learnModel_15 :
	mkdir ChromHMM_hg19_15states_output; \

	${java} -jar -Xmx10G ${ChromHMM}/ChromHMM.jar LearnModel -p 20 binarized_beds ChromHMM_hg19_15states_output 15 hg19

# run ChromHMM LearnModel to call 18 chromatin states
learnModel_18 :
	mkdir ChromHMM_hg19_18states_output2; \

	${java} -jar -Xmx10G ${ChromHMM}/ChromHMM.jar LearnModel -p 20 binarized_beds ChromHMM_hg19_18states_output2 18 hg19

chromhmm_18 :	binarizeBeds learnModel_18

# ------------------------------------------------------------------------------
#  ChromHMM on binarized FindER Peaks (fold > 2)
# ------------------------------------------------------------------------------
# generate 200bp binned bed files for each chromosome
chr_names = $(shell cat chr_names.tsv)

bin_200_loop = $(foreach chr, ${chr_names}, ChromHMM_fpeaks_dir/bin_200/$(chr).bed)
bin_200 :${bin_200_loop}

ChromHMM_fpeaks_dir/bin_200/%.bed :
	cat /projects/sdlee_prj/installs/ChromHMM/CHROMSIZES/hg19.txt | \
	awk '{ if($$1=="$*") print $$1"\t"0"\t"$$2 }' | awk '{ gsub("chr", ""); print }' | \
	${bedtools}/windowMaker -b - -w 200 > tmp/tmp.bed; \
	# remove final bin since it will be <200bp
	sed '$$d' tmp/tmp.bed > $@


is_bin_200_loop = $(foreach lib, ${libs_no_input}, ChromHMM_fpeaks_dir/$(lib)_is)
is_bin_200 : ${is_bin_200_loop}

ChromHMM_fpeaks_dir/%_is :
	mkdir ChromHMM_fpeaks_dir2; \
	mkdir ChromHMM_fpeaks_dir2/$*; \
	for chr in ${chr_names}; do \
		${bedtools}/intersectBed -wao -a ChromHMM_fpeaks_dir/bin_200/$$chr.bed \
		-b $*/peaks/peaks.sorted.bed > ChromHMM_fpeaks_dir2/$*/$$chr.bed; \
	done

intersect_bin_200_loop = $(foreach lib, ${libs_no_input}, ChromHMM_fpeaks_dir/$(lib)_intersect)
intersect_bin_200 : ${intersect_bin_200_loop}

ChromHMM_fpeaks_dir/%_intersect : 
	mkdir ChromHMM_fpeaks_dir/$*; \
	for chr in ${chr_names}; do \
		${bedtools}/intersectBed -wao -a ChromHMM_fpeaks_dir/bin_200/$$chr.bed \
		-b $*/peaks/peaks.fold2.bed > ChromHMM_fpeaks_dir/$*/$$chr.bed; \
	done

cell_lines ?= NHA NHAA2 NHAH9 F8 F8A2 F8E10

# generate binarized bed files from R first
# this loop will add the expected sample\tchromosome on the first line of each binarized file

add_first_line_loop = $(foreach line, ${cell_lines}, ChromHMM_fpeaks_dir/binarized_beds/$(line).bed)
add_first_line : ${add_first_line_loop}

ChromHMM_fpeaks_dir/binarized_beds/%.bed :
	for chr in ${chr_names}; do \
		echo -e "$*\t$${chr}\n$$(cat ChromHMM_fpeaks_dir2/binarized_beds/$*_$${chr}_binary.txt)" \
		> ChromHMM_fpeaks_dir2/binarized_beds/$*_$${chr}_binary.txt; \
	done

learnModel_fpeaks_18 :
	mkdir ChromHMM_fpeaks2_18; \

	${java} -jar -Xmx100G ${ChromHMM}/ChromHMM.jar LearnModel -p 80 ChromHMM_fpeaks_dir/binarized_beds \
	ChromHMM_fpeaks2_18 18 hg19

learnModel_fpeaks_15 :
	mkdir ChromHMM_fpeaks_15; \

	${java} -jar -Xmx100G ${ChromHMM}/ChromHMM.jar LearnModel -p 80 ChromHMM_fpeaks_dir/binarized_beds \
	ChromHMM_fpeaks_15 15 hg19

# ------------------------------------------------------------------------------
#  jaccard indeces between chromhmm segmentations
# ------------------------------------------------------------------------------
./chromhmm_18_states.tsv :
	for i in {1..18}; do echo $$i | awk '{print "E"$$i}' >> $@; done

chromhmm_18_states = $(shell cat ./chromhmm_18_states.tsv)

chromhmm_18_jac_loop = $(foreach state, ${chromhmm_18_states}, ./chromhmm_18_jacs2/$(state)/jacs.tsv)
chromhmm_18_jac : ${chromhmm_18_jac_loop}

./chromhmm_18_jacs2/%/jacs.tsv :
	mkdir ./chromhmm_18_jacs2/$*; \
	for i in ./ChromHMM_hg19_18states_output2/$*/*; do \
		for j in ./ChromHMM_hg19_18states_output2/$*/*; do \
			a=$$(echo $$i"__"$$j | awk '{gsub("./ChromHMM_hg19_18states_output2/$*/", ""); print}' | \
			     awk '{gsub("_18_segments.bed", ""); print}'); \
			b=$$(${bedtools}/bedtools jaccard -a $$i -b $$j | tail -1 | awk '{print $$3}'); \
			echo $${a}	$${b} >> $@; \
		done; \
	done

# ------------------------------------------------------------------------------
#  Generate bigWig files using deepTools (bin size = 10), normalize to RPKM
# ------------------------------------------------------------------------------
bam_to_bw_loop = $(foreach lib, ${libs}, $(lib)/bam/gsc.sorted.rpkm.bw)
bam_to_bw : ${bam_to_bw_loop}

%/bam/gsc.sorted.rpkm.bw :
	${deeptools}/bamCoverage -p 80 \
	-b $*/bam/gsc.sorted.bam \
	--binSize 10 \
	--normalizeUsing RPKM \
	--ignoreForNormalization chrX \
	--extendReads \
	-o $@

# normalize bigwig to input (default is log2 ratio)
norm_bw_loop = $(foreach pair, ${lib_pairs}, $(pair)/input_norm.bw)
norm_bw : ${norm_bw_loop}

%/input_norm.bw :
	ip=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
        input=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	${deeptools}/bigwigCompare -p 20 \
	-b1 $${ip}/bam/gsc.sorted.bam.bw \
	-b2 $${input}/bam/gsc.sorted.bam.bw \
	-o $${ip}/bam/input_norm.bw -of "bigwig"

nodupe_to_bw_loop = $(foreach lib, ${libs}, $(lib)/bam/gsc.nodupe.rpkm.bw)
nodupe_to_bw : ${nodupe_to_bw_loop}

%/bam/gsc.nodupe.rpkm.bw :
	${deeptools}/bamCoverage -p 80 \
	-b $*/bam/gsc.nodupe.bam \
	--binSize 10 \
	--normalizeUsing RPKM \
	--ignoreForNormalization chrX \
	--extendReads \
	-o $@

nodupe_log2_bw_loop = $(foreach pair, ${lib_pairs}, $(pair)/gsc.nodupe.log2.bw)
nodupe_log2_bw : ${nodupe_log2_bw_loop}

%/gsc.nodupe.log2.bw :
	ip=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	input=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	${deeptools}/bigwigCompare -p 100 \
	-b1 $${ip}/bam/gsc.nodupe.rpkm.bw \
	-b2 $${input}/bam/gsc.nodupe.rpkm.bw \
	--binSize 10 \
	-o $${ip}/bam/gsc.nodupe.log2.bw \
	--outFileFormat "bigwig"

# ------------------------------------------------------------------------------
#  Merge replicates (mean RPKM between rep1 & 2 into new bigWig)
# ------------------------------------------------------------------------------
rep_pairs = $(shell cat rep_pairs.tsv | grep -v "input")

mean_rpkm_bw = $(foreach pair, ${rep_pairs}, merged_bw/$(pair).bw)
mean_rpkm_bw_loop : ${mean_rpkm_bw}

merged_bw/%.bw :
	sample=$$(echo $* | awk '{ gsub("_rep.*", ""); print }'); \
	mark=$$(echo $* | awk '{ gsub(".*_", ""); print }'); \
	p1=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	p2=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	${deeptools}/bigwigCompare -p 80 \
	-b1 $${p1}/bam/gsc.sorted.rpkm.bw \
	-b2 $${p2}/bam/gsc.sorted.rpkm.bw \
	--operation mean \
	--binSize 10 \
	--skipNAs \
	-o merged_bw/$${sample}_$${mark}.rpkm.bw

mean_nodupe_rpkm_bw = $(foreach pair, ${rep_pairs}, merged_bw/$(pair).nodupe.bw)
mean_nodupe_rpkm_bw_loop : ${mean_nodupe_rpkm_bw}

merged_bw/%.nodupe.bw :
	sample=$$(echo $* | awk '{ gsub("_rep.*", ""); print }'); \
	mark=$$(echo $* | awk '{ gsub(".*_", ""); print }'); \
	p1=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	p2=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	${deeptools}/bigwigCompare -p 80 \
	-b1 $${p1}/bam/gsc.nodupe.rpkm.bw \
	-b2 $${p2}/bam/gsc.nodupe.rpkm.bw \
	--operation mean \
	--binSize 10 \
	--skipNAs \
	-o merged_bw/rpkm/$${mark}/$${sample}.rpkm.bw

# merged bw between CIC-KO lines
NHA_KO_merge_bw_loop = $(foreach mod, ${hmods}, merged_bw/rpkm/$(mod)/NHA_KO.rpkm.bw)
NHA_KO_merge_bw : ${NHA_KO_merge_bw_loop}

merged_bw/rpkm/%/NHA_KO.rpkm.bw :
	${deeptools}/bigwigCompare -p 60 \
	-b1 merged_bw/rpkm/$*/NHAA2.rpkm.bw \
	-b2 merged_bw/rpkm/$*/NHAH9.rpkm.bw \
	--operation mean \
	--binSize 10 \
	--skipNAs \
	-o $@

F8_KO_merge_bw_loop = $(foreach mod, ${hmods}, merged_bw/rpkm/$(mod)/F8_KO.rpkm.bw)
F8_KO_merge_bw : ${F8_KO_merge_bw_loop}

merged_bw/rpkm/%/F8_KO.rpkm.bw :
	${deeptools}/bigwigCompare -p 60 \
	-b1 merged_bw/rpkm/$*/F8A2.rpkm.bw \
	-b2 merged_bw/rpkm/$*/F8E10.rpkm.bw \
	--operation mean \
	--binSize 10 \
	--skipNAs \
	-o $@


input_pairs = $(shell cat rep_pairs.tsv | grep  "input")

mean_nodupe_rpkm_input_bw = $(foreach ipair, ${input_pairs}, merged_bw/input/$(ipair).nodupe.bw)
mean_nodupe_rpkm_input_bw_loop : ${mean_nodupe_rpkm_input_bw}

merged_bw/input/%.nodupe.bw :
	sample=$$(echo $* | awk '{ gsub("_rep.*", ""); print }'); \
	p1=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	p2=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	${deeptools}/bigwigCompare -p 80 \
	-b1 $${p1}/bam/gsc.nodupe.rpkm.bw \
	-b2 $${p2}/bam/gsc.nodupe.rpkm.bw \
	--operation mean \
	--binSize 10 \
	--skipNAs \
	-o merged_bw/rpkm/input/$${sample}.rpkm.bw

mean_log2_input_bw = $(foreach pair, ${rep_pairs}, merged_bw/input_norm/$(pair).log2.bw)
mean_log2_input_loop : ${mean_log2_input_bw}

merged_bw/input_norm/%.log2.bw :
	sample=$$(echo $* | awk '{ gsub("_rep.*", ""); print }'); \
	mark=$$(echo $* | awk '{ gsub(".*_", ""); print }'); \
	p1=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	p2=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	${deeptools}/bigwigCompare -p 40 \
	-b1 $${p1}/bam/input_norm.bw \
	-b2 $${p2}/bam/input_norm.bw \
	--operation mean \
	--skipNAs \
	-o merged_bw/input_norm/$${sample}_$${mark}.log2.bw

nodupe_log2_input_bw_loop = $(foreach pair, ${rep_pairs}, merged_bw/input_norm/$(pair).nodupe.log2.bw)
nodupe_log2_input_bw : ${nodupe_log2_input_bw_loop}

merged_bw/input_norm/%.nodupe.log2.bw :
	sample=$$(echo $* | awk '{ gsub("_rep.*", ""); print }'); \
	mark=$$(echo $* | awk '{ gsub(".*_", ""); print }'); \
	p1=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	p2=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	${deeptools}/bigwigCompare -p 80 \
	-b1 $${p1}/bam/gsc.nodupe.log2.bw \
	-b2 $${p2}/bam/gsc.nodupe.log2.bw \
	--operation mean \
	--skipNAs \
	-o merged_bw/input_norm/$${mark}/$${sample}.nodupe.log2.bw


# ------------------------------------------------------------------------------
#  computeMatrix
# ------------------------------------------------------------------------------
compute_matrix = $(foreach ip, ${ips}, $(ip)/peaks/peak_matrix.gz)
compute_matrixes : ${compute_matrix}

%/peaks/peak_matrix.gz :
	${deeptools}/computeMatrix reference-point -p 20 \
	--referencePoint center \
	-b 10000 -a 10000 \
	-R $*/peaks/peaks.sorted.bed \
	-S $*/bam/gsc.sorted.bam.bw \
	-o $@

comp_mat_pos = $(foreach ip, ${ips}, $(ip)/peaks/pos_peak_matrix.gz)
comp_mat_pos_all : ${comp_mat_pos}

%/peaks/pos_peak_matrix.gz :
	${deeptools}/computeMatrix reference-point -p 20 \
	--referencePoint center \
	-b 10000 -a 10000 \
	-R $*/peaks/cov_across_pos_genes.bed \
	-S $*/bam/gsc.sorted.bam.bw \
	-o $@

groups ?= IDH1 NHA both

cmat_DER_enh_loop = $(foreach group, ${groups}, enhancers/DER_enh_act.$(group).gz)
cmat_DER_enh : ${cmat_DER_enh_loop}

enhancers/DER_enh_act.%.gz :
	bigwigs=$$(ls merged_bw/rpkm/*/*.rpkm.bw | grep "h3k27ac" | \
	grep "NHA.rpkm.bw\|F8.rpkm.bw\|KO"); \
	names=$$(echo F8_KO.h3k27ac F8.h3k27ac NHA_KO.h3k27ac NHA.h3k27ac); \
	${deeptools}/computeMatrix reference-point -p 60 \
	--referencePoint center \
	-b 2000 -a 2000 \
	-R enhancers/DER_enh_act.$*.bed \
	-S $${bigwigs} \
	--samplesLabel $${names} \
	-o $@

DER_enh_hm_loop = $(foreach group, ${groups}, enhancers/DER_enh_act.$(group).pdf)
DER_enh_hm : ${DER_enh_hm_loop}

enhancers/DER_enh_act.%.pdf :
	${deeptools}/plotHeatmap \
	-m enhancers/DER_enh_act.$*.gz \
	--colorList white,limegreen \
	--whatToShow 'heatmap and colorbar' \
	-out $@

# ------------------------------------------------------------------------------
#  plotHeatmap
# ------------------------------------------------------------------------------
plot_heatmap = $(foreach ip, ${ips}, $(ip)/peaks/peak_heatmap.png)
plot_heatmaps : ${plot_heatmap}

%/peaks/peak_heatmap.png :
	${deeptools}/plotHeatmap \
	-m $*/peaks/peak_matrix.gz \
	-out $@ \
	--missingDataColor 1 \
	--colorList 'white,black'

plot_pos_heatmap = $(foreach ip, ${ips}, $(ip)/peaks/pos_peak_heatmap.png)
plot_pos_heatmaps : ${plot_pos_heatmap}

%/peaks/pos_peak_heatmap.png :
	${deeptools}/plotHeatmap \
	-m $*/peaks/pos_peak_matrix.gz \
	-out $@ \
	--missingDataColor 1 \
	--colorList 'white,black'



# ------------------------------------------------------------------------------
#  call primed, active and poised enhancers in CIC-WT cells
# ------------------------------------------------------------------------------
parental ?= NHA F8

# 1) Consolidate putative enhancer regions as previously defined (see J.Topham's thesis or D.Pellacani's paper)
merge_h3k4me1_loop = $(foreach line, ${cell_lines}, enhancers/$(line).h3k4me1.bed)
merge_h3k4me1 : ${merge_h3k4me1_loop}

# combine h3k4me1 peaks within 450bp of each other
# filter out regions <600bp
enhancers/%.h3k4me1.bed :
	${bedtools}/mergeBed -i merged_peaks/rep_pairs/$*_h3k4me1.fold2.bed -d 450 | \
	awk '{if (($$3 - $$2) >= 600) print}' | sort -k1,1 -k2,2n > $@

# 2) filter out H3K4me1 union peaks within 2kb of a TSS, remove duplicated rows with sort -u
filt_tss_loop = $(foreach line, ${cell_lines}, enhancers/$(line).h3k4me1.filt.bed)
filt_tss : ${filt_tss_loop}

enhancers/%.h3k4me1.filt.bed :
	${bedtools}/closestBed -d -a enhancers/$*.h3k4me1.bed \
	-b ${pd}/tss.bed | awk '{ if ($$8 > 2000 ) print $$1"\t"$$2"\t"$$3 }' | \
	sort -u | sort -k1,1 -k2,2n > enhancers/$*.h3k4me1.filt.bed

# intersect putative enhancers with DER peaks (seperately for up and down)
enh_WT_DER_loop = $(foreach line, ${parental}, enhancers/$(line).DER_enh.bed)
enh_WT_DER : ${enh_WT_DER_loop}

enhancers/%.DER_enh.bed :
	grep "down" DER_peaks/$*_DER_peaks.bed | \
	${bedtools}/intersectBed -wo -a enhancers/$*.h3k4me1.filt.bed \
	-b - > enhancers/$*.DER_enh.down.bed; \
	grep "up" DER_peaks/$*_DER_peaks.bed | \
	${bedtools}/intersectBed -wo -a enhancers/$*.h3k4me1.filt.bed \
	-b - > enhancers/$*.DER_enh.up.bed

# bin 5kb padded DER enhancers and get rpkm
DER_enh_5kb_bin_loop = $(foreach line, ${parental}, enhancers/$(line).DER_enh.5kb.100bin.bed)
DER_enh_5kb_bin : ${DER_enh_5kb_bin_loop}

enhancers/%.DER_enh.5kb.100bin.bed :
	${bedtools}/windowMaker -b enhancers/$*.DER_enh.down.5kb.bed -w 100 -i srcwinnum | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > enhancers/$*.DER_enh.down.5kb.100bin.bed; \
	${bedtools}/windowMaker -b enhancers/$*.DER_enh.up.5kb.bed -w 100 -i srcwinnum | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > enhancers/$*.DER_enh.up.5kb.100bin.bed

DER_enh_5kb_bin_rpkm_loop = $(foreach line, ${parental}, enhancers/$(line).DER_enh.5kb.100bin.rpkm.bed)
DER_enh_5kb_bin_rpkm : ${DER_enh_5kb_bin_rpkm_loop}

enhancers/%.DER_enh.5kb.100bin.rpkm.bed :
	bigwigs=$$(ls $**_h3*/bam/gsc.nodupe.rpkm.bw | grep -v "me3"); \
	names=$$(echo $${bigwigs} | awk '{ gsub("/bam/gsc.nodupe.rpkm.bw", ""); print}'); \
	for direction in down up; do \
		${deeptools}/multiBigwigSummary BED-file -p 20 \
		--BED enhancers/$*.DER_enh.$$direction.5kb.100bin.bed \
		--bwfiles $${bigwigs} \
		--labels $${names} \
		--outRawCounts enhancers/tmp.bed \
		--outFileName enhancers/$*.DER_enh.$$direction.5kb.100bin.npz; \
		awk 'NR>1' enhancers/tmp.bed | sort -k1,1 -k2,2n \
		> enhancers/$*.DER_enh.$$direction.5kb.100bin.rpkm.bed; \
	done


## overlap DER peaks with all enhancer regions rather than call states

# 3) overlap enhancers with other marks
enh_state_loop = $(foreach line, ${parental}, enhancers/$(line).act.bed)
enh_state : ${enh_state_loop}

# import the output files into R to categorize enhancers
# note that $*.act.bed will have some regions marked by both H3K27ac and 27me3 -> filter out in R
# since H3K4me1 peaks will have different coordinates between NHA and F8, used the merged coordinates as the enhancer id

enhancers/%.act.bed :
	cat enhancers/NHA.h3k4me1.filt.bed enhancers/F8.h3k4me1.filt.bed | sort -k1,1 -k2,2n | \
        ${bedtools}/mergeBed > enhancers/NHA_F8.h3k4me1.filt.bed; \
	${bedtools}/intersectBed -wao -a enhancers/NHA_F8.h3k4me1.filt.bed \
	-b enhancers/$*.h3k4me1.filt.bed | \
	awk '{ if ($$7 < 600) print $$4"\t"$$5"\t"$$6"\t"$$1":"$$2":"$$3 }' > enhancers/$*.unmarked.bed; \
	${bedtools}/intersectBed -wao -a enhancers/NHA_F8.h3k4me1.filt.bed \
	-b enhancers/$*.h3k4me1.filt.bed | \
	awk '{ if ($$7 >= 600) print $$4"\t"$$5"\t"$$6"\t"$$1":"$$2":"$$3 }' > enhancers/$*.primed.bed; \
	${bedtools}/intersectBed -wao -a enhancers/$*.primed.bed \
	-b merged_peaks/rep_pairs/$*_h3k27ac.fold2.bed | \
	awk '{ if ($$7 >= 300) print $$1"\t"$$2"\t"$$3"\t"$$4 }' > enhancers/$*.act.bed; \
	${bedtools}/intersectBed -wao -a enhancers/$*.primed.bed \
	-b merged_peaks/rep_pairs/$*_h3k27me3.fold2.bed | \
	awk '{ if ($$7 >= 300) print $$1"\t"$$2"\t"$$3"\t"$$4 }' > enhancers/$*.pois.bed; \
	${bedtools}/intersectBed -wao -a enhancers/$*.primed.bed \
	-b merged_peaks/rep_pairs/$*_h3k9me3.fold2.bed | \
	awk '{ if ($$7 >= 300) print $$1"\t"$$2"\t"$$3"\t"$$4 }' > enhancers/$*.het.bed

# 4) overlap DER peaks with annotated enhancers
#enh_WT_DER_loop = $(foreach line, ${parental}, enhancers/$(line).DER_peaks.bed)
#enh_WT_DER : ${enh_WT_DER_loop}

#enhancers/%.DER_peaks.bed :
#	${bedtools}/intersectBed -wao -a enhancers/pt_enh_states.bed \
#	-b DER_peaks/$*_DER_peaks.bed > $@

# bin 5kb padded enhancers
enh_states ?= off primed active poised

enh_5kb_bin_loop = $(foreach state, ${enh_states}, enhancers/$(state).5kb.100bin.bed)
enh_5kb_bin : ${enh_5kb_bin_loop}

enhancers/%.5kb.100bin.bed :
	${bedtools}/windowMaker -b enhancers/pt_enh.NHA_$*.5kb_pad.bed -w 100 -i srcwinnum | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > enhancers/pt_enh.NHA_$*.5kb_pad.100bp_bin.bed ; \
	${bedtools}/windowMaker -b enhancers/pt_enh.F8_$*.5kb_pad.bed -w 100 -i srcwinnum | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > enhancers/pt_enh.F8_$*.5kb_pad.100bp_bin.bed

enhancers/pt_enh.5kb_pad.100bp_bin.bed :
	${bedtools}/windowMaker -b enhancers/pt_enh.5kb_pad.bed -w 100 -i srcwinnum | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $@

enh_5kb_bin_spkm_NHA_loop = $(foreach state, ${enh_states}, enhancers/pt_enh.NHA_$(state).5kb_pad.100bp_bin.spkm.bed)
enh_5kb_bin_spkm_NHA : ${enh_5kb_bin_spkm_NHA_loop}

enhancers/pt_enh.NHA_%.5kb_pad.100bp_bin.spkm.bed :
	bigwigs=$$(ls */bam/gsc.nodupe.rpkm.bw | grep "NHA_" | grep "h3k4me1\|h3k27"); \
	names=$$(echo $${bigwigs} | awk '{ gsub("/bam/gsc.nodupe.rpkm.bw", ""); print }'); \
	${deeptools}/multiBigwigSummary BED-file -p 60 \
	--BED enhancers/pt_enh.NHA_$*.5kb_pad.100bp_bin.bed \
	--bwfiles $${bigwigs} \
	--labels $${names} \
	--outRawCounts enhancers/pt_enh.NHA_$*.5kb_pad.100bp_bin.spkm.tab \
	--outFileName enhancers/pt_enh.NHA_$*.5kb_pad.100bp_bin.spkm.npz; \
	awk 'NR>1' enhancers/pt_enh.NHA_$*.5kb_pad.100bp_bin.spkm.tab | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $@

enh_5kb_bin_spkm_F8_loop = $(foreach state, ${enh_states}, enhancers/pt_enh.F8_$(state).5kb_pad.100bp_bin.spkm.bed)
enh_5kb_bin_spkm_F8 : ${enh_5kb_bin_spkm_F8_loop}

enhancers/pt_enh.F8_%.5kb_pad.100bp_bin.spkm.bed :
	bigwigs=$$(ls */bam/gsc.nodupe.rpkm.bw | grep "F8_" | grep "h3k4me1\|h3k27"); \
	names=$$(echo $${bigwigs} | awk '{ gsub("/bam/gsc.nodupe.rpkm.bw", ""); print }'); \
	${deeptools}/multiBigwigSummary BED-file -p 60 \
	--BED enhancers/pt_enh.F8_$*.5kb_pad.100bp_bin.bed \
	--bwfiles $${bigwigs} \
	--labels $${names} \
	--outRawCounts enhancers/pt_enh.F8_$*.5kb_pad.100bp_bin.spkm.tab \
	--outFileName enhancers/pt_enh.F8_$*.5kb_pad.100bp_bin.spkm.npz; \
	awk 'NR>1' enhancers/pt_enh.F8_$*.5kb_pad.100bp_bin.spkm.tab | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $@

enh_5kb_bin_spkm : enh_5kb_bin_spkm_NHA enh_5kb_bin_spkm_F8

# ------------------------------------------------------------------------------
#  Do CIC peaks associate with DER enhancers?
# ------------------------------------------------------------------------------
enhancers/DER_enh_act_iCIC.bed :
	cat enhancers/DER_enh_act.bed | sort -k1,1 -k2,2n | \
	${bedtools}/closestBed -D "ref" -a - \
	-b ../tf_first_analysis/NHA_CIC_peaks_filt.bed > $@

# ------------------------------------------------------------------------------
#  merge h3k4me1 peaks
# ------------------------------------------------------------------------------
split_h3k4me1_loop = $(foreach line, ${cell_lines}, regions/$(line)/h3k4me1_split.bed)
split_h3k4me1 : ${split_h3k4me1_loop}

./regions/%/h3k4me1_split.bed :
	rm -R regions/$*; \
	mkdir regions/$*; \

	${bedtools}/intersectBed -wao -a regions/h3k4me1_merged.no_tss.bed \
	-b merged_peaks/rep_pairs/$*_h3k4me1.fpeaks.bed | \
	sort -k1,1 -k2,2n > regions/$*/h3k4me1_mapper.bed; \

	cat regions/$*/h3k4me1_mapper.bed | \
	awk '{ if($$4!=".") print $$4"\t"$$5"\t"$$6"\t"$$1":"$$2":"$$3 }' > $@

#  intersect h3k4me1 regions with enh marks
enh_marks ?= h3k4me1 h3k27ac h3k27me3 h3k9me3

h3k4me1_is_loop = $(foreach line, ${cell_lines}, regions/$(line)/h3k4me1_is.bed)
h3k4me1_is : ${h3k4me1_is_loop}

./regions/%/h3k4me1_is.bed :
	for mark in ${enh_marks}; do \
		${bedtools}/intersectBed -a merged_peaks/rep_pairs/$*_h3k4me1.fpeaks.bed \
		-b merged_peaks/rep_pairs/$*_$${mark}.fpeaks.bed | \
		${bedtools}/intersectBed -wao -a regions/h3k4me1_merged.no_tss.bed \
		-b - > regions/$*/h3k4me1_i_$${mark}.bed; \
	done

# get spkm in 5kb padded enhancers
enh_marks2 ?= h3k4me1 h3k27me3 h3k27ac

#comp_enh_5kb = $(foreach line, ${cell_lines}, regions/$(line)/h3k4me1_is.bed)
#h3k4me1_is : ${h3k4me1_is_loop}

regions/enh_5kb.h3k27me3.gz :
	${deeptools}/computeMatrix reference-point -p 80 \
	-S merged_bw/rpkm/h3k27me3/*.rpkm.bw \
	-R regions/enh_5kb.bed \
	--referencePoint center \
	-b 2500 -a 2500 \
	--binSize 100 \
	--smartLabels \
	--skipZeros \
	--missingDataAsZero \
	-o $@

regions/enh_5kb.h3k27me3.pdf :
	${deeptools}/plotHeatmap  \
	--matrixFile regions/enh_5kb.h3k27me3.gz \
	--outFileName $@ \
	--kmeans 4

regions/enh_5kb.h3k27ac.gz :
	${deeptools}/computeMatrix reference-point -p 80 \
	-S merged_bw/rpkm/h3k27ac/*.rpkm.bw \
	-R regions/enh_5kb.bed \
	--referencePoint center \
	-b 2500 -a 2500 \
	--binSize 100 \
	--smartLabels \
	--skipZeros \
	--missingDataAsZero \
	-o $@

regions/enh_5kb.h3k27ac.pdf :
	${deeptools}/plotHeatmap  \
	--matrixFile regions/enh_5kb.h3k27ac.gz \
	--outFileName $@ \
	--kmeans 4
	
# ------------------------------------------------------------------------------
#  Merge peaks across replicates
# ------------------------------------------------------------------------------
rep_pairs = $(shell cat rep_pairs.tsv | grep -v "input")

# Take the union of peaks between rep1 and 2
merge_rep_peaks_loop = $(foreach pair, ${rep_pairs}, merged_peaks/$(pair)_mq0.01.bed)
merge_rep_peaks : ${merge_rep_peaks_loop}

./merged_peaks/%_mq0.01.bed :
	p1=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	p2=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	cell_line=$$(echo $${p1} | awk '{ gsub("_.*_.*", ""); print }'); \
	mark=$$(echo $${p1} | awk '{ gsub(".*_rep.", ""); print }'); \
	${bedtools}/intersectBed -wo -a $${p1}/peaks/macs2_q0.01_peaks.narrowPeak \
	-b $${p2}/peaks/macs2_q0.01_peaks.narrowPeak > ./merged_peaks/$${cell_line}$${mark}_mq0.01.bed; \
	${bedops}/closest-features --dist --closest --delim '\t' \
	merged_peaks/$${cell_line}$${mark}_mq0.01.bed ${pd}/tss_pad_1500_500.bed > \
	merged_peaks/$${cell_line}$${mark}_mq0.01_tss.bed; \
	rm merged_peaks/$${cell_line}$${mark}_mq0.01.bed
	
# FindER
merge_rep_peaks_f_loop = $(foreach pair, ${rep_pairs}, merged_peaks/$(pair)_f.bed)
merge_rep_peaks_f : ${merge_rep_peaks_f_loop}

./merged_peaks/%_f.bed :
	p1=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	p2=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	cell_line=$$(echo $${p1} | awk '{ gsub("_.*_.*", ""); print }'); \
	mark=$$(echo $${p1} | awk '{ gsub(".*_rep.", ""); print }'); \
	cat $${p1}/peaks/peaks.sorted.bed $${p2}/peaks/peaks.sorted.bed | \
	sort -k1,1 -k2,2n | ${bedtools}/mergeBed > ./merged_peaks/$${cell_line}$${mark}_f.bed; \
	${bedops}/closest-features --dist --closest --delim '\t' \
	merged_peaks/$${cell_line}$${mark}_f.bed ${pd}/tss_pad_1500_500.bed > \
	merged_peaks/$${cell_line}$${mark}_f_tss.bed

## Obtain Histone mod Coverage (RPKM) in union peaks
# ------------------------------------------------------------------------------
#  computeMatrix
# ------------------------------------------------------------------------------
NHA_ips = $(shell cut -f1 ./libs.tsv | grep -v "^F8")
F8_ips = $(shell cut -f1 ./libs.tsv | grep -v "^NHA")

marks ?= h3k4me1 h3k4me3 h3k27ac h3k27me3 h3k36me3 h3k9me3 input
union_rpkm_loop = $(foreach mark, ${marks}, ./macs2_mpeaks/$(mark)/union_rpkm.npz)
union_rpkm : ${union_rpkm_loop}

# Obtain average rpkm within macs2 union peaks
./macs2_mpeaks/%/union_rpkm.npz :
	samples=$$(echo *$*/bam/gsc.sorted.rpkm.bw | awk '{ gsub("/bam/gsc.sorted.rpkm.bw", ""); print }'); \
	${deeptools}/multiBigwigSummary BED-file -p 40 \
	-b *$*/bam/gsc.sorted.rpkm.bw \
	--BED ./macs2_mpeaks/$*/union.bed \
	--labels $${samples} \
	--outRawCounts macs2_mpeaks/$*/union_rpkm.tab \
	-o $@

# Plot correlation heatmap using rpkm in union peaks
h_marks ?= h3k4me1 h3k4me3 h3k27ac h3k27me3 h3k36me3 h3k9me3
union_rpkm_cor_loop = $(foreach mark, ${h_marks}, ./macs2_mpeaks/$(mark)/union_cor.png)
union_rpkm_cor : ${union_rpkm_cor_loop}

./macs2_mpeaks/%/union_cor.png :
	${deeptools}/plotCorrelation \
	-in ./macs2_mpeaks/$*/union_rpkm.npz \
	--corMethod spearman \
	--skipZeros \
	--plotTitle "Spearman Correlation of RPKM in Peaks" \
	--whatToPlot heatmap \
	--colorMap RdYlBu \
	--plotNumbers \
	-o $@

macs2_mpeaks/h3k27me3/NHA_union_mat.gz :
	${deeptools}/computeMatrix scale-regions -p 40 \
	-S NHA*h3k27me3/bam/gsc.sorted.rpkm.bw \
	-R macs2_mpeaks/h3k27me3/union.bed \
	-m 10000 \
	--samplesLabel NHA*h3k27me3 \
	--skipZeros \
	-o $@ \
	--outFileSortedRegions macs2_mpeaks/h3k27me3/NHA_union.sorted.bed

macs2_mpeaks/h3k27me3/F8_union_mat.gz :
	${deeptools}/computeMatrix scale-regions -p 40 \
	-S F8*h3k27me3/bam/gsc.sorted.rpkm.bw \
	-R macs2_mpeaks/h3k27me3/union.bed \
	-m 10000 \
	--samplesLabel F8*h3k27me3 \
	--skipZeros \
	-o $@ \
	--outFileSortedRegions macs2_mpeaks/h3k27me3/F8_union.sorted.bed

macs2_mpeaks/h3k27me3/NHA_union.png :
	${deeptools}/plotHeatmap -m macs2_mpeaks/h3k27me3/NHA_union_mat.gz \
	--colorList 'black,white' \
	--perGroup \
	-out $@

macs2_mpeaks/h3k27me3/F8_union.png :
	${deeptools}/plotHeatmap -m macs2_mpeaks/h3k27me3/F8_union_mat.gz \
	--colorList 'black,white' \
	--perGroup \
	-out $@

#using input normalized bigwigs
macs2_mpeaks/h3k27me3/NHA_union_log2ratio_mat.gz :
	${deeptools}/computeMatrix scale-regions -p 40 \
	-S NHA*h3k27me3/bam/input_norm.bw \
	-R macs2_mpeaks/h3k27me3/union.bed \
	-m 10000 \
	--samplesLabel NHA*h3k27me3 \
	--skipZeros \
	-o $@

macs2_mpeaks/h3k27me3/NHA_union_log2ratio.png :
	${deeptools}/plotProfile -m macs2_mpeaks/h3k27me3/NHA_union_log2ratio_mat.gz \
	--perGroup \
	-out $@

macs2_mpeaks/h3k27me3/F8_union_log2ratio_mat.gz :
	${deeptools}/computeMatrix scale-regions -p 40 \
	-S F8*h3k27me3/bam/input_norm.bw \
	-R macs2_mpeaks/h3k27me3/union.bed \
	-m 10000 \
	--samplesLabel F8*h3k27me3 \
	--skipZeros \
	-o $@

macs2_mpeaks/h3k27me3/F8_union_log2ratio.png :
	${deeptools}/plotProfile -m macs2_mpeaks/h3k27me3/F8_union_log2ratio_mat.gz \
	--perGroup \
	-out $@

# ------------------------------------------------------------------------------
#  Obtain SPKM in 4kb bins genome wide and at promoters (2kb within ensembl TSS)
# ------------------------------------------------------------------------------
gw_4kb_spkm_loop = $(foreach mark, ${marks}, spkm/$(mark)/hg19_mappable.4kb_bin.spkm.tab)
gw_4kb_spkm : ${gw_4kb_spkm_loop}

# Obtain spkm within 4kb bins genome wide
spkm/%/hg19_mappable.4kb_bin.spkm.tab :
	mkdir spkm/$*; \
	samples=$$(echo *$*/bam/gsc.nodupe.rpkm.bw | awk '{ gsub("/bam/gsc.nodupe.rpkm.bw", ""); print }'); \
	${deeptools}/multiBigwigSummary BED-file -p 60 \
	-b *$*/bam/gsc.nodupe.rpkm.bw \
	--BED ${pd}/hg19_mappable.4kb_bin.bed \
	--labels $${samples} \
	--outRawCounts $@ \
	-o spkm/$*/hg19_mappable.4kb_bin.spkm.npz

promoter_spkm_loop = $(foreach mark, ${marks}, spkm/$(mark)/hg19_ensembl_promoter.spkm.tab)
promoter_spkm : ${promoter_spkm_loop}

# Obtain spkm within promoters
spkm/%/hg19_ensembl_promoter.spkm.tab :
	samples=$$(echo *$*/bam/gsc.nodupe.rpkm.bw | awk '{ gsub("/bam/gsc.nodupe.rpkm.bw", ""); print }'); \
	${deeptools}/multiBigwigSummary BED-file -p 60 \
	-b *$*/bam/gsc.nodupe.rpkm.bw \
	--BED ${pd}/hg19_ensembl_promoter.bed \
	--labels $${samples} \
	--outRawCounts $@ \
	-o spkm/$*/hg19_ensembl_promoter.spkm.npz

## Histone mod Profiles Across ens75 TSS
# ------------------------------------------------------------------------------
#  computeMatrix
# ------------------------------------------------------------------------------
NHA_ips = $(shell cut -f1 ./libs.tsv | grep -v "^F8")
F8_ips = $(shell cut -f1 ./libs.tsv | grep -v "^NHA")

NHA_DE_comp_mat_loop = $(foreach ip, ${NHA_ips}, $(ip)/peaks/DE_matrix.gz)
NHA_DE_comp_mat : ${NHA_DE_comp_mat_loop}

#%/peaks/DE_matrix.gz :

test.gz :
	${deeptools}/computeMatrix reference-point -p 20 \
	-S F8_rep1_h3k27me3/bam/input_norm.bw F8A2_rep1_h3k27me3/bam/input_norm.bw F8E10_rep1_h3k27me3/bam/input_norm.bw \
	-R ${pd}tss.bed \
	-b 5000 -a 5000 \
	--skipZeros \
	-o ./test.gz \
	--outFileSortedRegions ./test.sorted.bed

test.heatmap :
	${deeptools}/plotHeatmap -m test.gz \
	--colorList 'black,white' \
	--perGroup \
	-out test.png \

# ------------------------------------------------------------------------------
#  Histone mod profiles across jeme enhancers (NHA)
# ------------------------------------------------------------------------------
./objects/jeme_enhancers_NHA_padbin5kb.bed :
	cat ${pd}/jeme_enhancers_NHA.bed | \
	sort -k1,1 -k2,2n | \
	awk '{print $$1"\t"int(($$2 + $$3)/2)}' | \
	awk '{print $$1"\t"($$2 - 5000)"\t"($$2 + 5000)}' | \
	awk '$$2 > 0' | \
	${bedtools}/windowMaker -b - -w 100 -i srcwinnum | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $@

./objects/jeme_enhancers_NHA_padbin5kb_rpkm.npz :
	bigwigs=$$(ls */bam/gsc.sorted.rpkm.bw | grep -v input); \
	names=$$(echo $${bigwigs} | awk '{ gsub("/bam/gsc.sorted.rpkm.bw", ""); print }'); \
	${deeptools}/multiBigwigSummary BED-file -p 20 \
	--bwfiles $${bigwigs} \
	--BED ./objects/jeme_enhancers_NHA_padbin5kb.bed \
	--labels $${names} \
	--outRawCounts ./objects/jeme_enhancers_NHA_padbin5kb_rpkm.tab \
	-out $@

# ------------------------------------------------------------------------------
#  Plot heatmap for DE histone mod peaks
# ------------------------------------------------------------------------------
./F8_DE_h3k27me3.gz :
	${deeptools}/computeMatrix scale-regions -p 40 \
	-S  merged_bw/F8*_h3k27me3.rpkm.bw \
	-R F8_DE_h3k27me3.bed \
	-m 10000 \
	--samplesLabel F8A2 F8E10 F8 \
	--skipZeros \
	-o $@

F8_DE_h3k27me3.png :
	${deeptools}/plotHeatmap -m F8_DE_h3k27me3.gz \
	--sortRegions no \
	--colorList 'white, black' \
	-out $@

# ------------------------------------------------------------------------------
#  LOCKs
# ------------------------------------------------------------------------------
LOCKs/h3k27me3/f_locks_0.1_cov.tab :
	${deeptools}/multiBamSummary BED-file -p 30 \
	-b *h3k27me3/bam/gsc.sorted.bam \
	--extendReads \
	--ignoreDuplicates \
	--BED LOCKs/h3k27me3/f_locks_0.1.bed \
	--labels $$(echo *h3k27me3) \
	-o LOCKs/h3k27me3/f_locks_cov.npz \
	--outRawCounts $@

LOCKs/h3k9me3/f_locks_0.1_cov.tab :
	${deeptools}/multiBamSummary BED-file -p 30 \
	-b *h3k9me3/bam/gsc.sorted.bam \
	--extendReads \
	--ignoreDuplicates \
	--BED LOCKs/h3k9me3/f_locks_0.1.bed \
	--labels $$(echo *h3k9me3) \
	-o LOCKs/h3k9me3/f_locks_cov.npz \
	--outRawCounts $@

LOCKs/h3k36me3/f_locks_0.1_cov.tab :
	${deeptools}/multiBamSummary BED-file -p 30 \
	-b *h3k36me3/bam/gsc.sorted.bam \
	--extendReads \
	--ignoreDuplicates \
	--BED LOCKs/h3k36me3/f_locks_0.1.bed \
	--labels $$(echo *h3k36me3) \
	-o LOCKs/h3k36me3/f_locks_cov.npz \
	--outRawCounts $@

broad_hmod ?= h3k9me3 h3k27me3 h3k36me3

m_locks_cov_loop = $(foreach hmod, ${broad_hmod}, ./LOCKs/${hmod}/m_locks_0.1_cov.tab)
m_locks_cov : ${m_locks_cov_loop}

./LOCKs/%/m_locks_0.1_cov.tab :
	${deeptools}/multiBamSummary BED-file -p 40 \
	-b *$*/bam/gsc.sorted.bam \
	--extendReads \
	--ignoreDuplicates \
	--BED LOCKs/$*/m_locks_0.1.bed \
	--labels $$(echo *$*) \
	-o LOCKs/$*/m_locks_cov.npz \
	--outRawCounts $@


# ------------------------------------------------------------------------------
#  Plot heatmap for DE LOCKs
# ------------------------------------------------------------------------------
./LOCKs/h3k27me3/F8.gz :
	${deeptools}/computeMatrix scale-regions -p 40 \
	-S  merged_bw/input_norm/F8*_h3k27me3.log2.bw \
	-R ./LOCKs/h3k27me3/F8_lost_heatmap.bed ./LOCKs/h3k27me3/F8_gained_heatmap.bed ./LOCKs/h3k27me3/F8_ns_heatmap.bed \
	-m 10000 \
	--samplesLabel F8A2 F8E10 F8 \
	--skipZeros \
	-o $@

./LOCKs/h3k27me3/F8.png :
	${deeptools}/plotHeatmap -m ./LOCKs/h3k27me3/F8.gz \
	--sortRegions no \
	--colorList 'white, black' \
	-out $@

F8_LOCK_heatmap : ./LOCKs/h3k27me3/F8.gz ./LOCKs/h3k27me3/F8.png


# ------------------------------------------------------------------------------
#  H3K9me3 and H3K27me3 coverage at pericentromeres
# ------------------------------------------------------------------------------
# padbin 100000kb flanking either side of each centromere

../shared_objects/hg19_pc_100kb_padbin.bed :
	cat ${pd}/hg19_pc_100kb.bed | \
	/projects/sdlee_prj/installs/bedtools2/bin/windowMaker -b - -w 100 -i srcwinnum | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $@

./pericentromeres/pc_100kb_padbin_h3k9me3_rpkm.tab :
	${deeptools}/multiBigwigSummary BED-file -p 30 \
	--bwfiles *h3k9me3/bam/gsc.sorted.rpkm.bw \
	--BED ../shared_objects/hg19_pc_100kb_padbin.bed \
	--labels *h3k9me3 \
	--outRawCounts $@ \
	-out pericentromeres/pc_100kb_padbin_h3k9me3_rpkm.npz; \
	cat $@ | sort -k1,1 -k2,2n > pericentromeres/pc_100kb_padbin_h3k9me3_rpkm.bed

./pericentromeres/pc_100kb_padbin_h3k27me3_rpkm.tab :
	${deeptools}/multiBigwigSummary BED-file -p 30 \
	--bwfiles *h3k27me3/bam/gsc.sorted.rpkm.bw \
	--BED ../shared_objects/hg19_pc_100kb_padbin.bed \
	--labels *h3k27me3 \
	--outRawCounts $@ \
	-out pericentromeres/pc_100kb_padbin_h3k27me3_rpkm.npz; \
	cat $@ | sort -k1,1 -k2,2n > pericentromeres/pc_100kb_padbin_h3k27me3_rpkm.bed


# ------------------------------------------------------------------------------
#  RNA-seq Coverage within LOCKs
# ------------------------------------------------------------------------------
./LOCKs/h3k27me3/f_locks_0.1_rna_rpkm.tab :
	bigwigs=$$(ls ../rna_seq/merged_bams/*.bw); \
	names=$$(echo $${bigwigs} | \
		 awk '{ gsub("../rna_seq/merged_bams/", ""); print }' | \
		 awk '{ gsub(".bw", ""); print }'); \
	${deeptools}/multiBigwigSummary BED-file -p 60 \
	--bwfiles $${bigwigs} \
	--BED ./LOCKs/h3k27me3/f_locks_0.1.bed \
	--labels $${names} \
	--outRawCounts $@ \
	-out ./LOCKs/h3k27me3/f_locks_0.1.npz

# ------------------------------------------------------------------------------
#  Assess overlap between DE LOCKs and genomic regions (TADs, LADs)
# ------------------------------------------------------------------------------

# Generate control bed file with same size distribution	as DE LOCKs
# sort after make
./LOCKs/h3k27me3/NHA_DE_q0.01_ctrl.bed :
	${bedtools2}/shuffleBed -i ./LOCKs/h3k27me3/NHA_DE_q0.01.bed \
	-g ${pd}/hg19_chromsizes.txt \
	-excl ${pd}/consensusBlacklist.bed \
        -noOverlapping \
	-seed 12072013 \
	> $@

./LOCKs/h3k27me3/F8_DE_q0.01_ctrl.bed :
	${bedtools2}/shuffleBed -i ./LOCKs/h3k27me3/F8_DE_q0.01.bed \
	-g ${pd}/hg19_chromsizes.txt \
	-excl ${pd}/consensusBlacklist.bed \
	-noOverlapping \
	-seed 12072014 \
	> $@

# Generate control bed file with same size distribution as NHA TADs
NHA_TAD_ctrl :
	${bedtools2}/shuffleBed -i ${pd}/GSE105194_ENCFF306YQN_topologically_associated_domains_hg19.bed \
	-g ${pd}/hg19_chromsizes.txt \
	-excl ${pd}/consensusBlacklist.bed \
	-noOverlapping \
	-seed 1 \
	> ${pd}NHA_TAD_ctrl.bed

NHA_TAD_ctrl2 :
	${bedtools2}/shuffleBed -i ${pd}/NHA_TAD.bed \
	-g ${pd}/hg19_chromsizes.txt \
	-excl ${pd}/consensusBlacklist.bed \
	-noOverlapping \
	-seed 1 \
	> ${pd}NHA_TAD_ctrl2.bed

# Generate control bed file with same size distribution as LADs
LAD_ctrl :
	${bedtools2}/shuffleBed -i ${pd}/Fibroblast_LAD_hg19.bed \
	-g ${pd}/hg19_chromsizes.txt \
	-excl ${pd}/consensusBlacklist.bed \
	-noOverlapping \
	-seed 2 \
	> ${pd}LAD_ctrl.bed

# intersect LOCKs with TADs
./LOCKs/h3k27me3/DE_q0.01_i_TAD.bed :
	${bedtools2}/intersectBed -wo -a ./LOCKs/h3k27me3/NHA_DE_q0.01.bed \
	-b ${pd}/GSE105194_ENCFF306YQN_topologically_associated_domains_hg19.bed \
	> ./LOCKs/h3k27me3/NHA_DE_q0.01_i_TAD.bed; \

	${bedtools2}/intersectBed -wo -a ./LOCKs/h3k27me3/NHA_DE_q0.01_ctrl.bed \
	-b ${pd}/GSE105194_ENCFF306YQN_topologically_associated_domains_hg19.bed \
	> ./LOCKs/h3k27me3/NHA_DE_q0.01_ctrl_i_TAD.bed; \
	
	${bedtools2}/intersectBed -wo -a ./LOCKs/h3k27me3/F8_DE_q0.01.bed \
	-b ${pd}/GSE105194_ENCFF306YQN_topologically_associated_domains_hg19.bed \
	> ./LOCKs/h3k27me3/F8_DE_q0.01_i_TAD.bed; \

	${bedtools2}/intersectBed -wo -a ./LOCKs/h3k27me3/F8_DE_q0.01_ctrl.bed \
	-b ${pd}/GSE105194_ENCFF306YQN_topologically_associated_domains_hg19.bed \
	> ./LOCKs/h3k27me3/F8_DE_q0.01_ctrl_i_TAD.bed

# intersect LOCKs with LADs
LOCK_i_LADs :
	${bedtools2}/intersectBed -wo -a ./LOCKs/h3k27me3/NHA_DE_q0.01.bed \
	-b ${pd}/Fibroblast_LAD_hg19.bed > ./LOCKs/h3k27me3/NHA_DE_q0.01_i_LAD.bed; \

	${bedtools2}/intersectBed -wo -a ./LOCKs/h3k27me3/NHA_DE_q0.01.bed \
	-b ${pd}/LAD_ctrl.bed > ./LOCKs/h3k27me3/NHA_DE_q0.01_i_LAD_ctrl.bed; \
	
	${bedtools2}/intersectBed -wo -a ./LOCKs/h3k27me3/F8_DE_q0.01.bed \
	-b ${pd}/Fibroblast_LAD_hg19.bed > ./LOCKs/h3k27me3/F8_DE_q0.01_i_LAD.bed; \
	
	${bedtools2}/intersectBed -wo -a ./LOCKs/h3k27me3/F8_DE_q0.01.bed \
	-b ${pd}/LAD_ctrl.bed > ./LOCKs/h3k27me3/F8_DE_q0.01_i_LAD_ctrl.bed


parental ?= NHA F8

# jaccard index between LOCKs and TADs
./LOCKs/h3k27me3/TAD_jacs.tsv :
	${bedtools2}/bedtools jaccard -a ./LOCKs/h3k27me3/NHA_DE_q0.01.bed \
	-b ${pd}/NHA_TAD.bed | tail -1 >> $@; \

	${bedtools2}/bedtools jaccard -a ./LOCKs/h3k27me3/NHA_DE_q0.01.bed \
	-b ${pd}/NHA_TAD_ctrl.bed | tail -1 >> $@; \

	${bedtools2}/bedtools jaccard -a ./LOCKs/h3k27me3/F8_DE_q0.01.bed \
	-b ${pd}/NHA_TAD.bed | tail -1 >> $@; \

	${bedtools2}/bedtools jaccard -a ./LOCKs/h3k27me3/F8_DE_q0.01.bed \
	-b ${pd}/NHA_TAD_ctrl.bed | tail -1 >> $@

./LOCKs/h3k27me3/DE_TAD_jacs.tsv :
	for i in ${parental}; do \
		${bedtools2}/bedtools jaccard -a ./LOCKs/h3k27me3/$${i}_DE.bed \
		-b ${pd}/NHA_TAD.bed | tail -1 >> $@; \
		${bedtools2}/bedtools jaccard -a ./LOCKs/h3k27me3/$${i}_DE.bed \
		-b ${pd}/NHA_TAD_ctrl.bed | tail -1 >> $@; \
	done

./LOCKs/h3k27me3/DE_LAD_jacs.tsv :
	for i in ${parental}; do \
		${bedtools2}/bedtools jaccard -a ./LOCKs/h3k27me3/$${i}_DE.bed \
		-b ${pd}/Fibroblast_LAD_hg19.bed | tail -1 >> $@; \
		${bedtools2}/bedtools jaccard -a ./LOCKs/h3k27me3/$${i}_DE.bed \
		-b ${pd}/LAD_ctrl.bed | tail -1 >> $@; \
	done

# jaccard index between LOCKs and LADs
h3k27me3_LOCK_i_LAD_jac :
	${bedtools2}/bedtools jaccard -a ./LOCKs/h3k27me3/NHA_DE_q0.01.bed \
	-b ${pd}/Fibroblast_LAD_hg19.bed | tail -1 >> ./LOCKs/h3k27me3/LAD_jacs.tsv; \

	${bedtools2}/bedtools jaccard -a ./LOCKs/h3k27me3/NHA_DE_q0.01.bed \
	-b ${pd}/LAD_ctrl.bed | tail -1 >> ./LOCKs/h3k27me3/LAD_jacs.tsv; \

	${bedtools2}/bedtools jaccard -a ./LOCKs/h3k27me3/F8_DE_q0.01.bed \
	-b ${pd}/Fibroblast_LAD_hg19.bed | tail -1 >> ./LOCKs/h3k27me3/LAD_jacs.tsv; \
	
	${bedtools2}/bedtools jaccard -a ./LOCKs/h3k27me3/F8_DE_q0.01.bed \
	-b ${pd}/LAD_ctrl.bed | tail -1 >> ./LOCKs/h3k27me3/LAD_jacs.tsv

# ------------------------------------------------------------------------------
#  Promoter Analysis
# ------------------------------------------------------------------------------
# Promoters are defined as regions within 2kb of TSS of ensembl genes
# A histone mark needed to cover at least 400bp of the promoter to be considered

p_libs ?= $(shell cut -f1 libs.tsv | grep 'h3k27me3\|h3k4me3')

promoter_is_peak_loop = $(foreach lib, ${p_libs}, promoter2/$(lib).bed)
promoter_is_peak : ${promoter_is_peak_loop}

promoter2/%.bed :
	${bedtools}/intersectBed -wao -a ${pd}/hg19_EPD_promoter_2kb.bed \
	-b $*/peaks/peaks.sorted.bed > $@

# obtain intersect of h3k4me3 and h3k27me3 peaks for each sample and intersect with promoters
biv_is_loop = $(foreach sample, ${samples}, promoter2/$(sample)_biv.bed)
biv_is : ${biv_is_loop}

promoter2/%_biv.bed :
	${bedtools}/intersectBed -a $*_h3k4me3/peaks/peaks.sorted.bed \
	-b $*_h3k27me3/peaks/peaks.sorted.bed | \
	${bedtools}/intersectBed -wao -a ${pd}/hg19_EPD_promoter_2kb.bed \
	-b - > $@

promoter_is_all : promoter_is_peak biv_is

# obtain rpkm for h3k4me3 and h3k27me3 within promoters
promoter2/h3k4me3.rpkm.tab :
	names=$$(echo *h3k4me3 | awk '{ gsub("_h3k4me3", ""); print }'); \
	${deeptools}/multiBigwigSummary BED-file -p 80 \
	--BED ${pd}hg19_EPD_promoter_2kb.bed \
	--bwfiles *h3k4me3/bam/gsc.nodupe.rpkm.bw \
	--outRawCounts $@ \
	--labels $${names} \
	--outFileName promoter2/h3k4me3.rpkm.npz; \

promoter2/h3k27me3.rpkm.tab :
	names=$$(echo *h3k27me3 | awk '{ gsub("_h3k27me3", ""); print }'); \
	${deeptools}/multiBigwigSummary BED-file -p 80 \
	--BED ${pd}hg19_EPD_promoter_2kb.bed \
	--bwfiles *h3k27me3/bam/gsc.nodupe.rpkm.bw \
	--outRawCounts $@ \
	--labels $${names} \
	--outFileName promoter2/h3k27me3.rpkm.npz

rpkm_prom : promoter2/h3k4me3.rpkm.tab promoter2/h3k27me3.rpkm.tab

promoter2/h3k4me3.mrep.rpkm.tab :
	names=$$(echo merged_bw/rpkm/h3k4me3/*.rpkm.bw | awk '{ gsub("merged_bw/rpkm/h3k4me3/", ""); print }' | \
	awk '{ gsub(".rpkm.bw", ""); print }'); \
	${deeptools}/multiBigwigSummary BED-file -p 80 \
	--BED ${pd}hg19_EPD_promoter_2kb.bed \
	--bwfiles merged_bw/rpkm/h3k4me3/*.rpkm.bw \
	--outRawCounts $@ \
	--labels $${names} \
	--outFileName promoter2/h3k4me3.mrep.rpkm.npz; \

promoter2/h3k27me3.mrep.rpkm.tab :
	names=$$(echo merged_bw/rpkm/h3k27me3/*.rpkm.bw | awk '{ gsub("merged_bw/rpkm/h3k27me3/", ""); print }' | \
	awk '{ gsub(".rpkm.bw", ""); print }'); \
	${deeptools}/multiBigwigSummary BED-file -p 80 \
	--BED ${pd}hg19_EPD_promoter_2kb.bed \
	--bwfiles merged_bw/rpkm/h3k27me3/*.rpkm.bw \
	--outRawCounts $@ \
	--labels $${names} \
	--outFileName promoter2/h3k27me3.mrep.rpkm.npz

rpkm_mrep_prom : promoter2/h3k4me3.mrep.rpkm.tab promoter2/h3k27me3.mrep.rpkm.tab

# bin promoter state beds by 100bp
bin_prom_state_loop = $(foreach line, ${cell_lines}, promoter2/binned.$(line).bed)
bin_prom_state : ${bin_prom_state_loop}

promoter2/binned.%.bed :
	cat promoter2/h3k4me3/$*.bed | \
	${bedtools}/windowMaker -b - -w 100 -i srcwinnum > promoter2/h3k4me3/$*.5kb.bin.bed; \

	cat promoter2/h3k27me3/$*.bed | \
	${bedtools}/windowMaker -b - -w 100 -i srcwinnum > promoter2/h3k27me3/$*.5kb.bin.bed; \

	cat promoter2/bivalent/$*.bed | \
	${bedtools}/windowMaker -b - -w 100 -i srcwinnum > promoter2/bivalent/$*.5kb.bin.bed

# obtain spkm in binned promoters to generate metaplots in R
prom_states ?= h3k4me3 h3k27me3 bivalent

spkm_prom_state_loop = $(foreach line, ${cell_lines}, promoter2/$(line).spkm.bed)
spkm_prom_state : ${spkm_prom_state_loop}

promoter2/%.spkm.bed :
	for state in ${prom_states}; do \
		bigwigs=$$(ls merged_bw/rpkm/*/$*.rpkm.bw); \
		names=$$(echo $${bigwigs} | awk '{ gsub("merged_bw/rpkm/", ""); print }' | \
		awk '{ gsub("/$*.rpkm.bw", ""); print }'); \
		${deeptools}/multiBigwigSummary BED-file -p 80 \
		--BED promoter2/$${state}/$*.5kb.bin.bed \
		--bwfiles $${bigwigs} \
		--labels $${names} \
		--outRawCounts promoter2/$${state}/$*.5kb.bin.spkm.bed \
		--outFileName promoter2/$${state}/$*.5kb.bin.spkm.npz; \
		awk 'NR>1' promoter2/$${state}/$*.5kb.bin.spkm.bed | sort -k1,1 -k2,2n \
		> promoter2/$${state}/$*.5kb.bin.spkm.sorted.bed; \
	done
a :
	${deeptools}/multiBigwigSummary BED-file -p 80 \
	--BED promoter2/h3k27me3/$*.5kb.bin.bed \
	--bwfiles $${bigwigs} \
	--labels $${names} \
	--outRawCounts promoter2/h3k27me3/$*.5kb.bin.spkm.bed \
	--outFileName promoter2/h3k27me3/$*.5kb.bin.spkm.npz; \
	awk 'NR>1' promoter2/h3k27me3/$*.5kb.bin.spkm.bed | sort -k1,1 -k2,2n \
	> promoter2/h3k27me3/$*.5kb.bin.spkm.sorted.bed; \
	${deeptools}/multiBigwigSummary BED-file -p 80 \
	--BED promoter2/bivalent/$*.5kb.bin.bed \
	--bwfiles $${bigwigs} \
	--labels $${names} \
	--outRawCounts promoter2/bivalent/$*.5kb.bin.spkm.bed \
	--outFileName promoter2/bivalent/$*.5kb.bin.spkm.npz; \
	awk 'NR>1' promoter2/bivalent/$*.5kb.bin.spkm.bed | sort -k1,1 -k2,2n \
	> promoter2/bivalent/$*.5kb.bin.spkm.sorted.bed

# obtain log2 input signal in binned promoters to generate metplots in R
log2_prom_state_loop = $(foreach line, ${cell_lines}, promoter2/$(line).log2.bed)
log2_prom_state : ${log2_prom_state_loop}

promoter2/%.log2.bed :
	${deeptools}/multiBigwigSummary BED-file -p 80 \
	--BED promoter2/h3k4me3/$*.bin.bed \
	--bwfiles merged_bw/input_norm/h3k4me3/$*.nodupe.log2.bw merged_bw/input_norm/h3k27me3/$*.nodupe.log2.bw \
	--outRawCounts promoter2/h3k4me3/$*.bin.log2.bed \
	--outFileName promoter2/h3k4me3/$*.bin.log2.npz; \
	awk 'NR>1' promoter2/h3k4me3/$*.bin.log2.bed | sort -k1,1 -k2,2n \
	> promoter2/h3k4me3/$*.bin.log2.sorted.bed

	${deeptools}/multiBigwigSummary BED-file -p 80 \
	--BED promoter2/h3k27me3/$*.bin.bed \
	--bwfiles merged_bw/input_norm/h3k4me3/$*.nodupe.log2.bw merged_bw/input_norm/h3k27me3/$*.nodupe.log2.bw \
	--outRawCounts promoter2/h3k27me3/$*.bin.log2.bed \
	--outFileName promoter2/h3k27me3/$*.bin.log2.npz; \
	awk 'NR>1' promoter2/h3k27me3/$*.bin.log2.bed | sort -k1,1 -k2,2n \
	> promoter2/h3k27me3/$*.bin.log2.sorted.bed

	${deeptools}/multiBigwigSummary BED-file -p 80 \
	--BED promoter2/bivalent/$*.bin.bed \
	--bwfiles merged_bw/input_norm/h3k4me3/$*.nodupe.log2.bw merged_bw/input_norm/h3k27me3/$*.nodupe.log2.bw \
	--outRawCounts promoter2/bivalent/$*.bin.log2.bed \
	--outFileName promoter2/bivalent/$*.bin.log2.npz; \
	awk 'NR>1' promoter2/bivalent/$*.bin.log2.bed | sort -k1,1 -k2,2n \
	> promoter2/bivalent/$*.bin.log2.sorted.bed

# ------------------------------------------------------------------------------
#  Motif enrichment analysis on DER peaks assoc. w/ enhancers
# ------------------------------------------------------------------------------
enh_act_loss_motifs :
	for size in 500 1000; do \
		cat enhancers/NHA_DER_enh_act_loss.bed | \
		/projects/jtopham_prj/installs/homer/bin/findMotifsGenome.pl - \
		${pd}/hg19.fa motifs/homer/NHA_DER_enh_act_loss_$${size} \
		-p 40 -S 15 -size $${size} -mask; \
	done

enh_act_gain_motifs :
	for size in 500 1000; do \
		cat enhancers/NHA_DER_enh_act_gain.bed | \
		/projects/jtopham_prj/installs/homer/bin/findMotifsGenome.pl - \
		${pd}/hg19.fa motifs/homer/NHA_DER_enh_act_gain_$${size} \
		-p 40 -S 15 -size $${size} -mask; \
	done

enh_act_loss_DE_motifs :
	for size in 500 1000; do \
		cat enhancers/NHA_DER_enh_act_loss_DE.bed | \
		/projects/jtopham_prj/installs/homer/bin/findMotifsGenome.pl - \
		${pd}/hg19.fa motifs/homer/NHA_DER_enh_act_loss_DE_$${size} \
		-p 40 -S 15 -size $${size} -mask; \
	done

enh_act_gain_DE_motifs :
	for size in 500 1000; do \
		cat enhancers/NHA_DER_enh_act_gain_DE.bed | \
		/projects/jtopham_prj/installs/homer/bin/findMotifsGenome.pl - \
		${pd}/hg19.fa motifs/homer/NHA_DER_enh_act_gain_DE_$${size} \
		-p 40 -S 15 -size $${size} -mask; \
	done

enhancers/motifs/NHA_DER_enh_act :
	cat enhancers/NHA_DER_enh_act.bed | \
	/projects/jtopham_prj/installs/homer/bin/findMotifsGenome.pl - \
	${pd}/hg19.fa enhancers/motifs/NHA_DER_enh_act_500 -p 40 -S 15 -size 500 -mask

enhancers/motifs/F8_DER_enh_act :
	cat enhancers/F8_DER_enh_act.bed | \
	/projects/jtopham_prj/installs/homer/bin/findMotifsGenome.pl - \
	${pd}/hg19.fa enhancers/motifs/F8_DER_enh_act_500 -p 40 -S 15 -size 500 -mask

enhancers/motifs/IDH_DER_enh_act :
	cat enhancers/IDH_DER_enh_act.bed | \
	/projects/jtopham_prj/installs/homer/bin/findMotifsGenome.pl - \
	${pd}/hg19.fa enhancers/motifs/IDH_DER_enh_act_500 -p 40 -S 15 -size 500 -mask


enhancers/motifs/NHA_DER_enh_primed :
	cat enhancers/NHA_DER_enh_primed.bed | \
	/projects/jtopham_prj/installs/homer/bin/findMotifsGenome.pl - \
	${pd}/hg19.fa enhancers/motifs/NHA_DER_enh_primed_500 -p 40 -S 15 -size 500 -mask

enhancers/motifs/F8_DER_enh_primed :
	cat enhancers/F8_DER_enh_primed.bed | \
	/projects/jtopham_prj/installs/homer/bin/findMotifsGenome.pl - \
	${pd}/hg19.fa enhancers/motifs/F8_DER_enh_primed_500 -p 40 -S 15 -size 500 -mask

enhancers/motifs/IDH_DER_enh_primed :
	cat enhancers/IDH_DER_enh_primed.bed | \
	/projects/jtopham_prj/installs/homer/bin/findMotifsGenome.pl - \
	${pd}/hg19.fa enhancers/motifs/IDH_DER_enh_primed_500 -p 40 -S 15 -size 500 -mask

# motif analysis on enhancers specifically DER in NHA
enhancers/motifs/NHA_loss_DER_enh :
	cat enhancers/NHA_loss_DER_enh.bed | \
	/projects/jtopham_prj/installs/homer/bin/findMotifsGenome.pl - \
	${pd}/hg19.fa enhancers/motifs/NHA_loss_DER_enh_500 -p 40 -S 15 -size 500 -mask

enhancers/motifs/NHA_gain_DER_enh :
	cat enhancers/NHA_gain_DER_enh.bed | \
	/projects/jtopham_prj/installs/homer/bin/findMotifsGenome.pl - \
	${pd}/hg19.fa enhancers/motifs/NHA_gain_DER_enh_500 -p 40 -S 15 -size 500 -mask

enhancers/motifs/F8_loss_DER_enh :
	cat enhancers/F8_loss_DER_enh.bed | \
	/projects/jtopham_prj/installs/homer/bin/findMotifsGenome.pl - \
	${pd}/hg19.fa enhancers/motifs/F8_loss_DER_enh_500 -p 40 -S 15 -size 500 -mask

enhancers/motifs/F8_gain_DER_enh :
	cat enhancers/F8_gain_DER_enh.bed | \
	/projects/jtopham_prj/installs/homer/bin/findMotifsGenome.pl - \
	${pd}/hg19.fa enhancers/motifs/F8_gain_DER_enh_500 -p 40 -S 15 -size 500 -mask


# motif analysis on enhancers that exhibit CIC and IDH1 interactivity
interact_cats ?= concordant NHA_only F8_only_discordant

enh_act_interact_motif_loop = $(foreach cat, ${interact_cats}, enhancers/motifs/$(cat)_500)
enh_act_interact_motif : ${enh_act_interact_motif_loop}

enhancers/motifs/%_500 :
	cat enhancers/motifs/act_$*.bed | \
	/projects/jtopham_prj/installs/homer/bin/findMotifsGenome.pl - \
	${pd}/hg19.fa enhancers/motifs/act_$*_500 -p 40 -S 15 -size 500 -mask

motif_h3k27ac_loop = $(foreach set, ${subsets}, motifs/DE_h3k27ac/$(set)_knownResults.html)
motif_h3k27ac : ${motif_h3k27ac_loop}

motifs/DE_h3k27ac/%_knownResults.html :
	for size in 500 1000; do \
		mkdir motif/DE_h3k27ac/F8/$*/$${size}; \
		cat motif/DE_h3k27ac/F8/$*/q0.05_fc1.bed | \
		/projects/jtopham_prj/installs/homer/bin/findMotifsGenome.pl - \
		/projects/marralab_cic_prj/shared_objects/hg19.fa motif/DE_h3k27ac/F8/$*/$${size} \
		-p 10 -S 10 -N 20000 -size $${size} -len 5,6,7,8,9,10; \
	done

# ------------------------------------------------------------------------------
#  identify overlapping hmod peaks with CIC (high confidence) peaks
# ------------------------------------------------------------------------------
CIC_hmod_loop = $(foreach mark, ${hmods}, CIC_binding/CIC_i$(mark).bed)
CIC_hmod : ${CIC_hmod_loop}

CIC_binding/CIC_i%.bed :
	for line in NHA NHAA2 NHAH9 F8 F8A2 F8E10; do \
		${bedtools}/intersectBed -wao -a ../tf_first_analysis/NHA_CIC_hc_peaks.sorted.bed \
		-b merged_peaks/rep_pairs/$${line}_$*.fold2.bed > CIC_binding/$*/$${line}_CIC_hmod.bed; \
	done

CIC_binding/CIC_DER_peaks.bed :
	for line in NHA F8; do \
		${bedtools}/intersectBed -wao -a ../tf_first_analysis/NHA_CIC_hc_peaks.sorted.bed \
		-b DER_peaks/$${line}_DER_peaks.sorted.bed > CIC_binding/CIC_$${line}_DER_peaks.bed; \
	done
