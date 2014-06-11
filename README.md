######################################################################
# CancerPipeline 

This is a simple pipeline to automate the processing of cancer samples downstream of variant calling. 

The steps are as follows:

**Annotation**

- Filter variants according to coverage/allelfreq cutoffs.
- Annotate variants with Annovar and SNPeff.
- Transform sample vcfs into mafs.
- Unite all samples into a common output maf/vcf.

**MuSiC**

- Prepare input files for MuSiC (i.e. intersection of samples available as bam and samples in the united maf).
- Run all MuSiC tools sequentially.

**MutSigCV**

- Run MutSigCV on the united maf.

##Input

File locations and settings for the pipeline *mypipeline* are saved as a Python dictionary in *pipeline_config.py*.

##Output

For each unfiltered input *mysample*, the following outputs are created in a new folder named *./mypipeline/*:

- *mysample.vcf*: all variants that passed the filtering step
- *mysample.maf*: all filtered variants, in maf format
- *mysample_annotated.tsv*: all filtered variants with annotation from the vcf itself, Annovar and SNPeff

The following files are created for all unfiltered inputs together:

- *all_samples_mypipeline.maf*: all filtered variants in maf format 
- *all_samples_mypipeline.tsv*: all filtered and annotated variants

The folder *analysis_mypipeline* contains intermediate and output files for MuSiC and MutSigCV.
The list of significant genes created by is in *./analysis_mypipeline/output/mypipeline_mutsigcv.sig_genes.txt*

##Usage

> *"python CancerPipeline.py my_pipeline".*

##Configuration

The name and settings for each pipeline are configured in **CancerPipelineConfig.py**.

This file contains a Python dictionary for each pipeline.

**Location settings:**

Use absolute paths!

- root: the root folder of this pipeline 
- raw_vcf_folder: the folder containing the unfiltered vcfs
- bam_folder: the folder containing the bam and bam.bai files
- bed: the bed file for this panel  
- ref: path to the reference sequence (in fasta format)
- whitelist: optionally, supply a list of samples from the raw_vcf_folder to be processed 

**Filtering settings (optional)**
- min_cov: Minimum accepted coverage for a variant position
- min_varfreq: Minimum variant frequency for a variant position
- min_qual: Minimum variant frequency for a variant position 

**Run flags (optional)**
- vcf_type: select the type of input vcf (iontorrent/illumina_strelka) (default: iontorrent)
- cpus: number of CPUs to use in parallel (default: 1)
- functional_analysis: if False, don't run MutSigCV and MuSiC (default: True)
- version_numbers_not_in_blacklist: legacy flag, don't use 

##Dependencies

CoverageCheck expects a Linux system with Python 2.7 installed.
The Python packages numpy (1.8.1+) and pandas (0.13.1+) are expected as well. 

The pipeline depends on the installation of the following 3rd-party tools:

- **Annovar** (including the databases for all annotation sources, see Annovar documentation at http://www.openbioinformatics.org/annovar/annovar_download.html)
- **SNPeff** (v3.5)
- **vcf2maf** (github.com/ckandoth/vcf2maf)
- **MutSigCV** (v1.4)
- **Genome Music** (v0.4)

The install directories of these tools need to be set in ToolConfig.py.




