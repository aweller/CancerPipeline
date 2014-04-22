cancer_pipeline
===============

This is a simple pipeline to automate the processing of cancer samples downstream of variant calling. 
File locations and settings are saved as a Python dictionary in pipeline_config.py.

The steps are as follows:

- Annotation

Filter variants according to coverage/allelfreq cutoffs.
Annotate variants with Annovar and SNPeff.
Transform sample vcfs into mafs.
Unite all samples into a common output maf/vcf.

- MuSiC

Prepare input files for MuSiC (i.e. intersection of samples available as bam and samples in the united maf).
Run all MuSiC tools sequentially.

- MutSigCV

Run MutSigCV on the united maf.
