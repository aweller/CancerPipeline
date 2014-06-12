######################################################################
# CancerPipeline 

This is a simple pipeline to automate the processing of cancer samples downstream of variant calling. 

The pipeline steps are as follows:

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

CancerPipeline is written in Python/ruffus and aware of which steps have been run before on which files. Each rerun of the pipeline will thus only 
touch the files that have **changed** since the last time it was run, not stupidly rerun all of them.

##Input

The pipeline expects all unfiltered input vcf to be located in one folder.

    .
    └── root
        └── vcfs_raw
            ├── sample1.vcf
            └── sample2.vcf

Per pipeline run, it will then create an output folder for the filtered and annotated vcfs and an output folder for the functional analysis tools.

    .
    └── root
        ├── analysis_mypipeline
        │   ├── input
        │   │   └── mypipeline_input.txt
        │   └── output
        │       └── mypipeline_sign_genes.txt
        ├── vcfs_mypipeline
        │   ├── sample1_annotated.tsv
        │   ├── sample1.vcf
        │   ├── sample2_annotated.tsv
        │   └── sample2.vcf
        └── vcfs_raw
            ├── sample1.vcf
            └── sample2.vcf

File locations and settings per pipeline are saved into a regular text file in Python configuration format (see below).
The file contains information on each pipeline that was run, serving both as a config file for the current runs and a log file for past runs.

**Location settings**

Use absolute paths!

- root: the root folder of this pipeline 
- raw_vcf_folder: the folder containing the unfiltered vcfs
- bed: the bed file for this panel  
- ref: path to the reference sequence (in fasta format)
- bam_folder (optional): the folder containing the bam and bam.bai files. If not supplied, MuSiC will not be started.
- whitelist (optional): a list of samples from the raw_vcf_folder to process (default: use whole folder) 
- blacklist (optional): a list of samples from the raw_vcf_folder to NOT process (default: use whole folder) 

**Filtering settings (optional)**
- min_cov: Minimum accepted coverage for a variant position
- min_varfreq: Minimum variant frequency for a variant position
- min_qual: Minimum variant frequency for a variant position 

**Run flags (optional)**
- vcf_type: select the type of input vcf (iontorrent/illumina_strelka) (default: iontorrent)
- cpus: number of CPUs to use in parallel (default: 1)
- functional_analysis: if False, don't run MutSigCV and MuSiC (default: True)
- verbose_logging: if set to True, will result in more output while running (default: False)
- version_numbers_not_in_blacklist: legacy flag, don't use

**Example definition file**

    [Quasar_A1_VC1_VF03]
    
    root = /home/user/root/
    raw_vcf_folder = /home/user/root/vcfs_raw/
    bam_folder = /home/user/test/bams/
    bed = /home/user/root/data/mypipeline.bed
    ref = /home/user/root/data/hg19.fasta

    cpus = 1
    verbose_logging = False 
    min_cov = 100
    min_varfreq = 0.05


##Usage

If there's only one pipeline in the definition file, the only argument needed is the name of the file:

    > python CancerPipeline.py pipeline_definitions.txt

If there's more than one pipeline defined, we need to select a pipeline with the 2nd argument::

    > python CancerPipeline.py pipeline_definitions.txt mypipeline

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

##Dependencies

CoverageCheck was developed on Ubuntu 13.10 with Python 2.7 and the Python packages numpy (1.8.1+), pandas (0.13.1+) and ruffus (2.4+). 

The pipeline depends on the installation of the following 3rd-party tools:

- **Annovar** (including the databases for all annotation sources, see Annovar documentation at http://www.openbioinformatics.org/annovar/annovar_download.html)
- **SNPeff** (v3.5)
- **vcf2maf** (github.com/ckandoth/vcf2maf)
- **MutSigCV** (v1.4)
- **Genome Music** (v0.4)

The install directories of these tools need to be set in **ToolConfig.py**.




