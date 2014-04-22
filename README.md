cancer_pipeline
===============

The ruffus pipeline from a bunch of vcfs to annotation and downstream analysis in MuSiC and MutSigCV.

This is a simple pipeline to automate the processing of cancer samples downstream of variant calling. 

The steps are as follows:

-------------------------
Filter variants according to coverage/allelfreq cutoffs.
Annotate variants with Annovar and SNPeff.
Transform sample vcfs into mafs.
Unite all samples into a common output maf/vcf.
-------------------------
Prepare input files for MuSiC (i.e. intersection of samples available as bam and samples in the united maf).
Run all MuSiC tools sequentially.
-------------------------
Run MutSigCV on the united maf.
