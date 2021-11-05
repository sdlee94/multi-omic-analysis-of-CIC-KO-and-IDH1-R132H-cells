wd = $(shell pwd)
pd = /projects/marralab_cic_prj/shared_objects/
PERCENT := %

bedtools = /projects/sdlee_prj/installs/bedtools2/bin
bedops = /projects/jtopham_prj/installs/bedops/bin
deeptools = /gsc/software/linux-x86_64-centos6/python-2.7.12/bin
defiant = /projects/sdlee_prj/installs/defiant/defiant
finder = /home/mbilenky/bin/Solexa_Java/FindER.1.0.1e.jar
java = /gsc/software/linux-x86_64/jre1.8.0_66/bin/java
UCSC = /projects/sdlee_prj/installs/UCSC_tools

# ------------------------------------------------------------------------------
#  set up directory
# ------------------------------------------------------------------------------

libs = $(shell cut -f 1 libs.tsv)
lib_locs = $(shell cat libs.tsv | awk '{print$$1"__"$$2}')

setupdir = $(foreach lib, ${libs}, $(lib)/bam)
dir : ${setupdir}

# create directory structure
%/bam :
	mkdir $*; mkdir $*/bam

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
flagstats = $(foreach lib, ${libs}, $(lib)/bam/flagstat.txt)

sortbams : ${filter_sort}
flags : ${flagstats}

# sort and index
%/bam/gsc.sorted.bam :
	samtools sort $*/bam/gsc.bam $*/bam/gsc.sorted; \
	samtools index $*/bam/gsc.sorted.bam \
	$*/bam/gsc.sorted.bam.bai

# flagstat
%/bam/flagstat.txt :
	samtools flagstat $*/bam/gsc.bam > $@
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#  filter CpGs for those with >5 coverage and sort
# ------------------------------------------------------------------------------
sort_cpg_loop = $(foreach lib, ${libs}, $(lib)/novo5mC/CpG.c5.sorted.bed)
sort_cpg : $(sort_cpg_loop)

%/novo5mC/CpG.c5.sorted.bed :
	cat $*/novo5mC/gsc.q5.5mC.CpG.autoChrXY.sum_cpg.quant | \
	awk '{if($$4+$$5>5) print $$1"\t"$$2"\t"$$2+1"\t"$$4"\t"$$5"\t"$$7}' | \
	${bedops}/sort-bed --max-mem 100G --tmpdir ./tmp - > $@

# ------------------------------------------------------------------------------
# Filter CpGs for those with greater than 5 coverage
# ------------------------------------------------------------------------------
filter_cpg = $(foreach lib, ${libs}, $(lib)/novo5mC/gsc.q5.5mC.CpG.autoChrXY.sum_cpg.quant.c5)
filter_cpgs : $(filter_cpg)

%/novo5mC/gsc.q5.5mC.CpG.autoChrXY.sum_cpg.quant.c5 : %/novo5mC/gsc.q5.5mC.CpG.autoChrXY.sum_cpg.quant
	cat $< | awk '{if($$4+$$5>5) print $$0}' > $@

# ------------------------------------------------------------------------------
#  subset CpGs with >=5 coverage for CpG islands and shores
# ------------------------------------------------------------------------------
island_cpg_loop = $(foreach lib, ${libs}, $(lib)/novo5mC/CpG.c5.cpgi.bed)
island_cpg : $(island_cpg_loop)

%/novo5mC/CpG.c5.cpgi.bed :
	${bedtools}/intersectBed -a $*/novo5mC/CpG.c5.sorted.bed \
	-b ref_data/cpg_islands.sorted.bed > $@;
	${bedtools}/intersectBed -a $*/novo5mC/CpG.c5.sorted.bed \
	-b ref_data/cpg_shores.sorted.bed > $*/novo5mC/CpG.c5.cpgsh.bed

# ------------------------------------------------------------------------------
# Reformat methylation calls into appropriate input for methylKit (R)
# Generic format: chrBase        chr    base     strand  coverage  freqC  freqT
# 		  chr21.9764539  chr21  9764539  R       12        25.00  75.00
# ------------------------------------------------------------------------------
reformat_me = $(foreach lib, ${libs}, $(lib)/novo5mC/gsc.q5.5mC.CpG.autoChrXY.sum_cpg.quant.mekit)
reformat_all_me : $(reformat_me)
%/novo5mC/gsc.q5.5mC.CpG.autoChrXY.sum_cpg.quant.mekit : %/novo5mC/gsc.q5.5mC.CpG.autoChrXY.sum_cpg.quant
	less $< | awk '{sum=$$4+$$5; freq_C=$$5/sum*100; freq_T=$$4/sum*100; print "chr"$$1"."$$2"\tchr"$$1"\t"$$2"\t"$$3"\t"sum"\t"freq_C"\t"freq_T}' > $@

# ------------------------------------------------------------------------------
#  Find closest genes to differentially methylated regions
# ------------------------------------------------------------------------------

# First run methyl_diff.sh and DMR.dynamic.sh

#DMR_libs = $(shell cut -f 1 DM_libs.tsv)

#testa:
#	echo ${DMR_libs}

#bedtools_closest = $(foreach lib, ${DMR_libs}, $(lib)/DMR.s300.c3.gene.bed)
#find_closest_genes : ${bedtools_closest}

#%/DMR.s300.c3.gene.bed :
#	${bedtools}/closestBed -a $*/DMR..s300.c3 -b ${pd}/Homo_sapiens_GRCh37_genes.sorted.bed -d \
#	| awk '{print $$1"\t"$$2"\t"$$3"\t"$$5"\t"$$6"\t"$$7"\t"$$9"\t"$$10"\t"$$12"\t"$$13"\t"$$14"\t"$$15}' > $@

#test_loop : $(foreach lib, ${DMR_libs}, $(lib)/check.loop.txt)
#test_a : ${test_loop}

#%/check.loop.txt :
#	echo $* > $@

# ------------------------------------------------------------------------------
#  Intersect DMRs with CpG island (using files generated from R)
# ------------------------------------------------------------------------------
#wgbs_dmrs_CpG_islands.bed :
#	${bedtools}/intersectBed -a wgbs_dmrs.bed -b cpg_island_genes_df.bed > $@

#intersect_dm = $(foreach lib, ${DMR_libs}, $(lib)/wgbs_dm_CpGs_in_dm_islands.bed)
#intersect_dm_all : ${intersect_dm}

#wgbs_dm_CpGs_in_dm_islands.bed :
#	${bedtools}/intersectBed -wo -a wgbs_dm_cpgs.bed -b $*/wgbs_dmrs_CpG_islands.bed \
#	| awk '{print $$1"\t"$$2"\t"$$3"\t"$$5"\t"$$6"\t"$$7"\t"$$18"\t"$$20}' > $@

# ------------------------------------------------------------------------------
#  Average CpG methylation in genomic features
# ------------------------------------------------------------------------------
# CpG Islands
avg_cpg_in_cpgi_loop = $(foreach lib, ${libs}, ./genomic_features/CpG_islands/$(lib).mean_cpg.bed)
avg_cpg_in_cpgi : ${avg_cpg_in_cpgi_loop}

./genomic_features/CpG_islands/%.mean_cpg.bed :
	cat $*/novo5mC/cpg_c5.bed | sort -k1,1 -k2,2n > $*/novo5mC/cpg_c5.sorted.bed; \
	${bedtools2}/mapBed -a ./ref_data/cpg_islands.sorted.bed \
	-b $*/novo5mC/cpg_c5.sorted.bed -c 4 -o mean > $@

avg_cpg_in_cpg_shores_loop = $(foreach lib, ${libs}, ./genomic_features/CpG_shores/$(lib).mean_cpg.bed)
avg_cpg_in_cpg_shores : ${avg_cpg_in_cpg_shores_loop}

./genomic_features/CpG_shores/%.mean_cpg.bed :
	${bedtools2}/mapBed -a ./ref_data/cpg_shores.sorted.bed \
	-b $*/novo5mC/cpg_c5.sorted.bed -c 4 -o mean > $@

avg_cpg_at_tss_loop = $(foreach lib, ${libs}, ./genomic_features/tss_1500/$(lib).mean_cpg.bed)
avg_cpg_at_tss : ${avg_cpg_at_tss_loop}

./genomic_features/tss_1500/%.mean_cpg.bed :
	${bedtools2}/mapBed -a ./ref_data/refseq_tss_1500.sorted.bed \
	-b $*/novo5mC/cpg_c5.sorted.bed -c 4 -o mean > $@; \
	${bedtools2}/mapBed -a ./ref_data/refseq_tss_200.sorted.bed \
	-b $*/novo5mC/cpg_c5.sorted.bed -c 4 -o mean > ./genomic_features/tss_200/$*.mean_cpg.bed

avg_cpg_gene_loop = $(foreach lib, ${libs}, ./genomic_features/5_utr/$(lib).mean_cpg.bed)
avg_cpg_gene : ${avg_cpg_gene_loop}

./genomic_features/5_utr/%.mean_cpg.bed :
	${bedtools2}/mapBed -a ./ref_data/refseq_5_utr.sorted.bed \
	-b $*/novo5mC/cpg_c5.sorted.bed -c 4 -o mean > $@; \
	${bedtools2}/mapBed -a ./ref_data/refseq_gene_body.sorted.bed \
	-b $*/novo5mC/cpg_c5.sorted.bed -c 4 -o mean > ./genomic_features/gene_body/$*.mean_cpg.bed; \
	${bedtools2}/mapBed -a ./ref_data/refseq_3_utr.sorted.bed \
	-b $*/novo5mC/cpg_c5.sorted.bed -c 4 -o mean > ./genomic_features/3_utr/$*.mean_cpg.bed

avg_cpg_intergenic_loop = $(foreach lib, ${libs}, ./genomic_features/intergenic/$(lib).mean_cpg.bed)
avg_cpg_intergenic : ${avg_cpg_intergenic_loop}

./genomic_features/intergenic/%.mean_cpg.bed :
	${bedtools2}/mapBed -a ${pd}/hg19_intergenic.bed \
	-b $*/novo5mC/cpg_c5.sorted.bed -c 4 -o mean > $@

avg_cpg_CIC_target_tss_loop = $(foreach lib, ${libs}, ./genomic_features/CIC_target/$(lib).mean_cpg.bed)
avg_cpg_CIC_target_tss : ${avg_cpg_CIC_target_tss_loop}

./genomic_features/CIC_target/%.mean_cpg.bed :
	${bedtools2}/mapBed -a ../Weissman_et_al/Weiss_G144_CIC_targets_tss.sorted.bed \
	-b $*/novo5mC/cpg_c5.sorted.bed -c 4 -o mean > $@


# ------------------------------------------------------------------------------
#  Average CpG methylation across genomic bins
# ------------------------------------------------------------------------------
./objects/hg19_bin2kb.bed :
	cat objects/cpgs_range_hg19.bed | \
	${bedtools}/windowMaker -b - -w 2000 -i srcwinnum | sort -k1,1 -k2,2n > $@

./objects/hg19_bin2kb.bed2 :
	cat objects/cpgs_range_hg19.bed | \
	${bedtools}/windowMaker -b - -w 2000 -i srcwinnum | \
	${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - > $@

./objects/hg19_bin10kb.bed :
	cat objects/cpgs_range_hg19.bed | \
	${bedtools}/windowMaker -b - -w 10000 -i srcwinnum | \
	${bedops}/sort-bed --max-mem 80G --tmpdir ./tmp - > $@

intersect_cpg_10kb_loop = $(foreach lib, ${libs}, genome_10kb/$(lib)/cpg_c5_10kb.bed)
intersect_cpg_10kb : $(intersect_cpg_10kb_loop)

genome_10kb/%/cpg_c5_10kb.bed :
	mkdir genome_10kb/$*; \
	${bedtools}/intersectBed -wao -a objects/hg19_bin10kb.bed \
	-b $*/novo5mC/CpG.c5.sorted.bed > $@

avg_cpg_2kb_loop = $(foreach lib, ${libs}, $(lib)/novo5mC/avg_cpg_2kb_c5.bed)
avg_cpg_2kb : $(avg_cpg_2kb_loop)

%/novo5mC/avg_cpg_2kb_c5.bed :
	cat $*/novo5mC/gsc.q5.5mC.CpG.autoChrXY.sum_cpg.quant.c5 | \
	awk '{ print ($$1"\t"$$2"\t"$$2"\t"$$7) }' | \
	sort -k1,1 -k2,2n > tmp/$*_cpg_c5.bed; \
	${bedtools2}/mapBed -a objects/hg19_bin2kb.bed -b tmp/$*_cpg_c5.bed -c 4 -o mean > $@

objects/test_avg.bed :
	${bedtools}/mapBed -a objects/hg19_bin2kb.bed -b tmp/cpg_c5.bed -c 4 -o mean > $@

tmp/cpg_c5.bed :
	cat F8_rep1/novo5mC/gsc.q5.5mC.CpG.autoChrXY.sum_cpg.quant.c5 | \
	awk '{ print ($$1"\t"$$2"\t"$$2"\t"$$7) }' > $@

avg_cpg_in_bed_loop = $(foreach lib, ${libs}, ./LOCKs/h3k27me3/f_locks_0.1_$(lib).cpg_mean.bed)
avg_cpg_in_bed : ${avg_cpg_in_bed_loop}

./LOCKs/h3k27me3/f_locks_0.1_%.cpg_mean.bed :
	cat ../histone_mod_chip/LOCKs/h3k27me3/f_locks_0.1.bed | \
	sort -k1,1 -k2,2n > ./LOCKs/h3k27me3/f_locks_0.1.bed; \
	${bedtools2}/mapBed -a ./LOCKs/h3k27me3/f_locks_0.1.bed \
	-b tmp/$*_cpg_c5.bed -c 4 -o mean > $@

# ------------------------------------------------------------------------------
#  DMR Analysis using Defiant
# ------------------------------------------------------------------------------
# https://github.com/hhg7/defiant

run_defiant :
	mekit=$$(echo novo5mC/gsc.q5.5mC.CpG.autoChrXY.sum_cpg.quant.mekit); \
	${defiant} -c 10 -cpu 50  -L WT,KO -l NHAA2 -N -v \
	-i NHA_rep1/$${mekit},NHA_rep2/$${mekit} NHAA2_rep1/$${mekit},NHAA2_rep2/$${mekit}; \
	${defiant} -c 10 -cpu 50  -L WT,KO -l NHAH9 -N -v \
	-i NHA_rep1/$${mekit},NHA_rep2/$${mekit} NHAH9_rep1/$${mekit},NHAH9_rep2/$${mekit}; \
	${defiant} -c 10 -cpu 50  -L WT,KO -l F8A2 -N -v \
	-i F8_rep1/$${mekit},F8_rep2/$${mekit} F8A2_rep1/$${mekit},F8A2_rep2/$${mekit}; \
	${defiant} -c 10 -cpu 50  -L WT,KO -l F8E10 -N -v \
	-i F8_rep1/$${mekit},F8_rep2/$${mekit} F8E10_rep1/$${mekit},F8E10_rep2/$${mekit}; \


defiant_wt :
	mekit=$$(echo novo5mC/gsc.q5.5mC.CpG.autoChrXY.sum_cpg.quant.mekit); \
	${defiant} -c 10 -cpu 50  -L WT,KO -l IDH -N -v \
	-i NHA_rep1/$${mekit},NHA_rep2/$${mekit} F8_rep1/$${mekit},F8_rep2/$${mekit}; \


#intersect defiant results b/w CIC-KO replicate lines and merge DMRs
intersect_DMRs_NHA :
	${bedtools}/intersectBed -wo -a thesis_DMRs/NHAA2.bed \
	-b thesis_DMRs/NHAH9.bed > thesis_DMRs/NHA.i.bed; \
	${bedtools}/intersectBed -wo -a thesis_DMRs/F8A2.bed \
	-b thesis_DMRs/F8E10.bed > thesis_DMRs/F8.i.bed

merge_DMRs :
	for cell in NHA F8; do \
		cut -f-4,9 thesis_DMRs/concordant_$$cell.i.bed > tmp/p1; \
		cut -f5- thesis_DMRs/concordant_$$cell.i.bed > tmp/p2; \
		cat tmp/p1 tmp/p2 | sort -k1,1 -k2,2n | \
		${bedtools}/mergeBed -d 1 > thesis_DMRs/$$cell.merged.bed; \
	done

# get original DMR info from each CIC-KO cell line for each merged DMR
merged_DMR_info :
	for cell in NHA F8; do \
		${bedtools}/intersectBed -wo -a thesis_DMRs/$$cell.merged.bed \
		-b thesis_DMRs/concordant_$$cell.i.bed > thesis_DMRs/$$cell.merged.info.bed; \
	done

# get all CpGs from each cell line for each merged DMR
DMR_cpgs :
	for cell in NHA F8; do \
		${bedtools}/intersectBed -wo -a thesis_DMRs/$$cell.merged.bed \
		-b $$cell_rep*/novo5mC/CpG.c5.sorted.bed > thesis_DMRs/$$cell.merged.CpGs.bed; \
	done

# ------------------------------------------------------------------------------
#  Pad DMRs +/- 5kb from centre, bin by 200bp and obtain methylation values
# ------------------------------------------------------------------------------
# too few DMRs in F8

reference_lines ?= NHA IDH1

DMR_5kb_bin_loop = $(foreach line, ${reference_lines}, thesis_DMRs/$(line).5kb.200bin.bed)
DMR_5kb_bin : $(DMR_5kb_bin_loop)

thesis_DMRs/%.5kb.200bin.bed :
	for direction in hypo hyper; do \
		cat thesis_DMRs/$*.$$direction.5kb.bed | \
		${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - | \
		${bedtools}/windowMaker -b - -w 200 -i srcwinnum | \
		${bedops}/sort-bed --max-mem 20G --tmpdir ./tmp - \
		> thesis_DMRs/$*.$$direction.5kb.200bin.bed; \
	done

thesis_DMRs/NHA.5kb.200bin.fme.bed :
	for direction in hypo hyper; do \
		for cell in NHA NHAA2 NHAH9; do \
			${bedtools}/mapBed -a thesis_DMRs/NHA.$$direction.5kb.200bin.bed \
			-b merged_rep/$$cell.CpG.c5.bed -c 4 -o mean \
			> thesis_DMRs/NHA.$$direction.5kb.200bin.$$cell.fme.bed; \
		done; \
	done

thesis_DMRs/IDH1.5kb.200bin.fme.bed :
	for direction in hypo hyper; do \
		for cell in NHA F8; do \
			${bedtools}/mapBed -a thesis_DMRs/IDH1.$$direction.5kb.200bin.bed \
			-b merged_rep/$$cell.CpG.c5.bed -c 4 -o mean \
			> thesis_DMRs/IDH1.$$direction.5kb.200bin.$$cell.fme.bed; \
		done; \
	done

DMR_bins_fme :
	for cell in NHAA2 NHAH9 F8A2 F8E10; do \
		${bedtools}/mapBed -a thesis_DMRs/$$cell.5kb.200_bin.bed \
		-b merged_rep/$$cell.CpG.c5.bed -c 4 -o mean > thesis_DMRs/$$cell.5kb.200_bin.$$cell.fme.bed; \
	done; \
	for cell in NHAA2 NHAH9; do \
		${bedtools}/mapBed -a thesis_DMRs/$$cell.5kb.200_bin.bed \
		-b merged_rep/NHA.CpG.c5.bed -c 4 -o mean > thesis_DMRs/$$cell.5kb.200_bin.NHA.fme.bed; \
	done; \
	for cell in F8A2 F8E10; do \
		${bedtools}/mapBed -a thesis_DMRs/$$cell.5kb.200_bin.bed \
		-b merged_rep/F8.CpG.c5.bed -c 4 -o mean > thesis_DMRs/$$cell.5kb.200_bin.F8.fme.bed; \
	done; \
	for cell in NHA F8; do \
		${bedtools}/mapBed -a thesis_DMRs/IDH1.5kb.200_bin.bed \
		-b merged_rep/$$cell.CpG.c5.bed -c 4 -o mean > thesis_DMRs/IDH1.5kb.200_bin.$$cell.fme.bed; \
	done

# ------------------------------------------------------------------------------
#  prepare bed files for MethylSeekR
#  columns are as follows: chr	start	end	C+T	C
# ------------------------------------------------------------------------------
libs ?= NHA_rep1 NHA_rep2 NHAA2_rep1 NHAA2_rep2 NHAH9_rep1 NHAH9_rep2 \
F8_rep1 F8_rep2 F8A2_rep1 F8A2_rep2 F8E10_rep1 F8E10_rep2

MSR_bed_loop = $(foreach lib, ${libs}, MethylSeekR/$(lib)/cpg.bed)
MSR_bed : ${MSR_bed_loop}

MethylSeekR/%/cpg.bed :
	mkdir MethylSeekR/$*; \
	cat $*/novo5mC/gsc.q5.5mC.CpG.autoChrXY.sum_cpg.quant | \
	awk '{ print $$1"\t"$$2"\t"$$2"\t"$$4+$$5"\t"$$5 }' > $@; \
	sort -k1,1 -k2,2n $@ > MethylSeekR/$*/cpg.sorted.bed

# ------------------------------------------------------------------------------
#  generate bed files for each MethylSeekR segment type
# ------------------------------------------------------------------------------
cell_lines ?= NHA NHAA2 NHAH9 F8 F8A2 F8E10

FMR_loop = $(foreach line, ${cell_lines}, MethylSeekR/$(line)/FMR.bed)
FMR : ${FMR_loop}

MethylSeekR/%/FMR.bed :
	${bedtools}/subtractBed -a ../shared_objects/hg19_mappable2.bed \
	-b MethylSeekR/$*/PMD.bed > tmp/tmpPMD.bed; \
	${bedtools}/subtractBed -a tmp/tmpPMD.bed \
	-b MethylSeekR/$*/LMR.bed > tmp/tmpLMR.bed; \
	${bedtools}/subtractBed -a tmp/tmpLMR.bed \
	-b MethylSeekR/$*/UMR.bed > $@

# ------------------------------------------------------------------------------
#  intersect CpGs with MethylSeekR segment beds
# ------------------------------------------------------------------------------
intersect_cpg_MSR_loop = $(foreach lib, ${libs}, MethylSeekR/$(lib)/MSR.bed)
intersect_cpg_MSR : ${intersect_cpg_MSR_loop}

MethylSeekR/%/MSR.bed :
	line=$$(echo $* | awk '{ gsub("_rep.*", ""); print}'); \
	${bedtools}/intersectBed -wao -a MethylSeekR/$${line}/FMR.bed \
	-b $*/novo5mC/CpG.c5.sorted.bed > MethylSeekR/$*/cpg_FMR.bed; \
	${bedtools}/intersectBed -wao -a MethylSeekR/$${line}/UMR.bed \
	-b $*/novo5mC/CpG.c5.sorted.bed > MethylSeekR/$*/cpg_UMR.bed; \
	${bedtools}/intersectBed -wao -a MethylSeekR/$${line}/LMR.bed \
	-b $*/novo5mC/CpG.c5.sorted.bed > MethylSeekR/$*/cpg_LMR.bed; \
	${bedtools}/intersectBed -wao -a MethylSeekR/$${line}/PMD.bed \
	-b $*/novo5mC/CpG.c5.sorted.bed > MethylSeekR/$*/cpg_PMD.bed; \

# ------------------------------------------------------------------------------
#  Generate BigWig files for IGV Visualization
# ------------------------------------------------------------------------------
cpg_bg_loop = $(foreach lib, ${libs}, $(lib)/bam/cpg.bw)
cpg_bg : ${cpg_bg_loop}

%/bam/cpg.bw :
	cat $*/novo5mC/gsc.q5.5mC.CpG.autoChrXY.sum_cpg.quant | \
	awk '{ print "chr"$$1"\t"$$2"\t"$$2"\t"$$7*100 }' | sort -k1,1 -k2,2n > $*/novo5mC/cpg.bedgraph; \
	${UCSC}/bedGraphToBigWig $*/novo5mC/cpg.bedgraph ${pd}/hg19_canon_chromsizes.txt $@

# ------------------------------------------------------------------------------
#  Merge BigWigs between replicates
# ------------------------------------------------------------------------------
rep_pairs = $(shell cat rep_pairs.tsv)

cpg_bg_repmerge_loop = $(foreach pair, ${rep_pairs}, $(pair)/bam/cpg.repmerge.bw)
cpg_bg_repmerge : ${cpg_bg_repmerge_loop}

%/bam/cpg.repmerge.bw :
	name=$$(echo $* | awk '{ gsub("_rep.*", ""); print }'); \
	r1=$$(echo $* | awk '{ gsub("__.*", ""); print }'); \
	r2=$$(echo $* | awk '{ gsub(".*__", ""); print }'); \
	${deeptools}/bigwigCompare -p 30 \
	-b1 $${r1}/bam/cpg.bw \
	-b2 $${r2}/bam/cpg.bw \
	--binSize 1 \
	--operation mean \
	--skipNAs \
	-o merged_rep/$${name}.cpg.bw

# ------------------------------------------------------------------------------
#  intersect DMRs with genomic features
# ------------------------------------------------------------------------------
DMR_gfeat :
	for cell in NHA F8; do \
		${bedtools}/intersectBed -wao -a thesis_DMRs/$$cell.merged.bed \
		-b ${pd}/hg19_gfeatures.bed | uniq > thesis_DMRs/$$cell.merged.gfeat.bed; \
	done

IDH1_DMR_gfeat :
	cut -f-3 thesis_DMRs/IDH1.bed | \
	${bedtools}/intersectBed -wao -a - \
	-b ${pd}/hg19_gfeatures.bed > thesis_DMRs/IDH1.merged.gfeat.bed

# ------------------------------------------------------------------------------
#  assess overlap b/w DMRs
# ------------------------------------------------------------------------------
overlap_DMRs :
	${bedtools}/intersectBed -wo -a thesis_DMRs/NHA.merged.bed \
	-b thesis_DMRs/F8.merged.bed > NHA.F8.bed; \
	${bedtools}/intersectBed -wo -a thesis_DMRs/NHA.merged.bed \
	-b thesis_DMRs/IDH1.bed > NHA.IDH1.bed; \
	${bedtools}/intersectBed -wo -a thesis_DMRs/F8.merged.bed \
	-b thesis_DMRs/IDH1.bed > F8.IDH1.bed
 

# ------------------------------------------------------------------------------
#  assess overlap b/w DMRs and CIC binding sites
# ------------------------------------------------------------------------------
overap_DMR_CIC_peaks :
	${bedtools}/intersectBed -wo -a thesis_DMRs/NHA.merged.info.bed \
	-b ../tf_first_analysis/NHA_CIC_hc_peaks.bed > thesis_DMRs/NHA.CIC.bed; \
	${bedtools}/intersectBed -wo -a thesis_DMRs/F8.merged.info.bed \
	-b ../tf_first_analysis/NHA_CIC_hc_peaks.bed > thesis_DMRs/F8.CIC.bed


# ------------------------------------------------------------------------------
#  calculate avg CpGme at CIC binding sites
# ------------------------------------------------------------------------------
avg_cpg_CIC_peaks :
	for NHA_line in NHA NHAA2 NHAH9; do \
		${bedtools}/mapBed -a ../tf_first_analysis/NHA_CIC_hc_peaks.sorted.bed \
		-b merged_rep/$$NHA_line.CpG.c5.bed -c 4 -o mean > CIC_hc_peaks.$$NHA_line.avg_cpgme.bed; \
	done; \
	for F8_line in F8 F8A2 F8E10; do \
		${bedtools}/mapBed -a ../tf_first_analysis/NHA_CIC_hc_peaks.sorted.bed \
		-b merged_rep/$$F8_line.CpG.c5.bed -c 4 -o mean > CIC_hc_peaks.$$F8_line.avg_cpgme.bed; \
	done

# ------------------------------------------------------------------------------
#  calculate avg CpGme +/- 2kb around TSS of CIC targets
# ------------------------------------------------------------------------------
avg_cpg_CIC_target_TSS :
	for NHA_line in NHA NHAA2 NHAH9; do \
		${bedtools}/mapBed -a CIC_target_TSS.2kb_pad.bed \
		-b merged_rep/$$NHA_line.CpG.c5.bed -c 4 -o mean > CIC_target_TSS.$$NHA_line.avg_cpgme.bed; \
	done; \
	for F8_line in F8 F8A2 F8E10; do \
		${bedtools}/mapBed -a CIC_target_TSS.2kb_pad.bed \
		-b merged_rep/$$F8_line.CpG.c5.bed -c 4 -o mean > CIC_target_TSS.$$F8_line.avg_cpgme.bed; \
	done

