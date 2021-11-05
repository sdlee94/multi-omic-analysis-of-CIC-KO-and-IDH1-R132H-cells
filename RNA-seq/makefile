# ------------------------------------------------------------------------------
#                           _      _      _
#              		 __(.)< __(.)> __(.)=
#         	         \___)  \___)  \___)
# ------------------------------------------------------------------------------
# the purpose of this pipeline is to work with the CIC RNA-seq data.
#
# ------------------------------------------------------------------------------
#              	            _      _      _
#          	    	 __(.)< __(.)> __(.)=
#          	    	 \___)  \___)  \___)
# ------------------------------------------------------------------------------

wd = $(shell pwd)
pd = /projects/marralab_cic_prj/shared_objects
PERCENT := %

bedtools = /projects/sdlee_prj/installs/bedtools2/bin
deeptools = /gsc/software/linux-x86_64-centos6/python-2.7.12/bin
subread = /projects/jtopham_prj/installs/subread-1.4.6-p5-source/bin
gtf_file = /projects/marralab_cic_prj/shared_objects/Homo_sapiens.GRCh37.75.gtf
leafcutter = /projects/sdlee_prj/installs/leafcutter
merge_reps = /projects/vporter_prj/tools/other/merge_bigwig.sh

# ------------------------------------------------------------------------------
#  set up directory
# ------------------------------------------------------------------------------

libs = $(shell cat libs.tsv | grep "/" | cut -f 1)
lib_locs = $(shell cat libs.tsv | awk '{print$$1"__"$$2}')

setupdir = $(foreach lib, ${libs}, $(lib)/bam)
dir : ${setupdir}

# create directory structure
%/bam :
	mkdir $*; mkdir $*/bam;	mkdir $*/rna;

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

# filter unmapped reads; sort and index
%/bam/gsc.sorted.bam :
	samtools view -bF 4 $*/bam/gsc.bam > $*/bam/tmp.bam; \
	samtools sort $*/bam/tmp.bam $*/bam/gsc.sorted; \
	samtools index $*/bam/gsc.sorted.bam \
	$*/bam/gsc.sorted.bam.bai; rm $*/bam/tmp.bam

# flagstat
%/bam/flagstat.txt :
	samtools flagstat $*/bam/gsc.bam > $@

# filter non-canoncal chromsomes
filter_canon_loop = $(foreach lib, ${libs}, $(lib)/bam/gsc.sorted.canon.bam)
filter_canon : ${filter_canon_loop}

%/bam/gsc.sorted.canon.bam :
	samtools view -o $@ $*/bam/gsc.sorted.bam ${pd}/canon_chr.txt

# ------------------------------------------------------------------------------
#  generate feature counts using SUBREAD
# ------------------------------------------------------------------------------

run_counts = $(foreach lib, ${libs}, $(lib)/rna/counts.tsv)
counts : ${run_counts}

%/rna/counts.tsv :
	cd $*/rna; \
	${subread}/featureCounts -pC -T 20 -t 'gene' -a ${gtf_file} -o counts.tsv \
	${wd}/$*/bam/gsc.sorted.bam

# ------------------------------------------------------------------------------
#  Generate bigWig files using deepTools (bin size = 10), normalize to RPKM
# ------------------------------------------------------------------------------
bam_to_bw_loop = $(foreach lib, ${libs}, $(lib)/bam/gsc.sorted.rpkm.bw)
bam_to_bw : ${bam_to_bw_loop}

%/bam/gsc.sorted.rpkm.bw :
	${deeptools}/bamCoverage -p 40 \
	-b $*/bam/gsc.sorted.canon.bam \
	--binSize 10 \
	--normalizeUsing RPKM \
	--ignoreForNormalization chrX \
	--extendReads \
	-o $@

# ------------------------------------------------------------------------------
#  Merge replicate BAMs and generate bigWigs to get mean RPKM across replicates
# ------------------------------------------------------------------------------
cell_lines ?= NHA NHAA2 NHAH9 F8 F8A2 F8E10 L54 L54N4D6 L54N8F4

# merge replicate BAMs and generate BAM index files
merge_bam_loop = $(foreach line, ${cell_lines}, ./merged_bams/$(line).bam)
merge_bam : ${merge_bam_loop}

./merged_bams/%.bam :
	samtools merge $@ \
	$*_rep1/bam/gsc.sorted.bam $*_rep2/bam/gsc.sorted.bam $*_rep3/bam/gsc.sorted.bam;
	samtools index $@ ./merged_bams/$*.bai

mbam_to_bw_loop = $(foreach line, ${cell_lines}, ./merged_bams/$(line).bw)
mbam_to_bw : ${mbam_to_bw_loop}

./merged_bams/%.bw :
	${deeptools}/bamCoverage -p 40 \
	-b ./merged_bams/$*.bam \
	--binSize 10 \
	--normalizeUsing RPKM \
	--ignoreForNormalization chrX \
	--extendReads \
	-o $@

# ------------------------------------------------------------------------------
#  Merge replicates (mean RPKM between rep1 & 2  & 3 into new bigWig)
# ------------------------------------------------------------------------------
cell_lines ?= NHA NHAA2 NHAH9 F8 F8A2 F8E10 L54 L54N4D6 L54N8F4

mean_rpkm_bw_loop = $(foreach cell_line, ${cell_lines}, merged_bw/$(cell_line).rpkm.bw)
mean_rpkm_bw : ${mean_rpkm_bw_loop}

./merged_bw/%.rpkm.bw :
	${merge_reps} -g 64 $@ \
	${pd}/hg19_chrsizes.txt \
	$*_rep1/bam/gsc.sorted.rpkm.bw \
	$*_rep2/bam/gsc.sorted.rpkm.bw \
	$*_rep3/bam/gsc.sorted.rpkm.bw

# ------------------------------------------------------------------------------
#  motif enrichment analysis using HOMER
# ------------------------------------------------------------------------------

NHA_knownResults_q0.05_up.html  :
	/projects/jtopham_prj/installs/homer/bin/findMotifs.pl ./results/NHA_res_q0.05_up_genes.txt \
	human ./results/homer/$@ -start -500 -end 500 -p 10
