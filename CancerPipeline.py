from ruffus import *
import sys
import os
import AnnotationRowParsers as AP
import SampleAnnotation as SA
import automate_vcf2maf as vcf2maf
import prepare_and_run_music as music
import automate_mutsig as mutsig
import CancerPipelineConfig as config
import pprint
import shutil
import logging
import re

# TODO
# - fix the TODO related to not having "_v1" in the current bam names 
# 
# - there are still "non-confident" variants in the filtered output
# - some buggy columns in the annotated files (e.g. numbers popping up in Codon_Change)

#configure logging 
logging.basicConfig(level=logging.DEBUG,
                format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                datefmt='%d-%m-%y %H:%M',
                #filename= "pipeline_log.txt",
                stream= sys.stdout)
                #filemode='w')

logging.debug("Logging DEBUG activated.")

########################################################################################
# Parse the config file  ###############################################################
########################################################################################

project_name = sys.argv[1]

config = config.config

if not config.get(project_name):
    print "Unknown project, please choose from the following or define a new project:"
    print "\n".join( config.keys() )
    sys.exit()
    
config = config[project_name]
config["name"] = project_name

root = config["root"]
#config["bam_folder"] = root + "bams/"
config["vcf_folder"] = root + "vcfs_%s/" % (project_name)
config["analysis_folder"] = root + "analysis_%s/" % (project_name)

if not os.path.exists(config["vcf_folder"]):
    os.mkdir(config["vcf_folder"])

if not os.path.exists(config["analysis_folder"]):
    os.mkdir(config["analysis_folder"])
    os.mkdir(config["analysis_folder"] + "input")
    os.mkdir(config["analysis_folder"] + "output")
    os.mkdir(config["analysis_folder"] + "output/roi_covgs/")
    os.mkdir(config["analysis_folder"] + "output/gene_covgs/")

# copy existing covg files into the current analysis folder
for root, dirs, files in os.walk(root):
    for name in files:
        if name.endswith(".covg"):
            source = root+"/"+name
            dest = config["analysis_folder"] + "output/roi_covgs/"
            if dest not in source:
                shutil.copy(source, dest)   

###########################################################
# Clean up the folders of intermediate files that might be left over from previous pipeline runs


def clean_folders_of_intermediate_files():

    dirty_folders = [config["vcf_folder"], config["raw_vcf_folder"]]
    dirty_words = ["annoin", "hg19", "snpeff", "refGene"] # intermediate files from Annovar/SNPeff

    for folder in dirty_folders:
        for filename in os.listdir(folder):
            for word in dirty_words:
                if word in filename:
                    os.remove(folder + filename)


########################################################################################
# Prepare the unfiltered vcfs ##########################################################
########################################################################################

###########################################################
# Annotation

raw_vcfs = [config["raw_vcf_folder"]+x for x in os.listdir(config["raw_vcf_folder"])
                  if x.endswith(".vcf") and
                  "QC" not in x and
                  "all" not in x and
                  "snpeff" not in x and
                  "anno" not in x]

logging.debug("Found %s raw input files:" % len(raw_vcfs))
logging.debug(raw_vcfs[:5])

if not config.get("version_numbers_not_in_blacklist"):
    # The usual case: both blacklist and samples contains "_v1":
    if config.get("blacklist"):
        logging.debug("Using blacklist %s" % config.get("blacklist"))
        blacklist = [x.strip() for x in open(config.get("blacklist")).readlines()]
        raw_vcfs = [x for x in raw_vcfs if x.split("/")[-1] not in blacklist] 
    if config.get("whitelist"):
        logging.debug("Using whitelist %s" % config.get("whitelist"))
        whitelist = [x.strip() for x in open(config.get("whitelist")).readlines()]
        raw_vcfs = [x for x in raw_vcfs if x.split("/")[-1] in whitelist] 

else:
    #TODO: remove the sub once the bam names are fixed
    # Blacklist doesnt contain "_v1" while the filesnames do, so they need to be trimmed before comparison
    if config.get("blacklist"):
        blacklist = [x.strip() for x in open(config.get("blacklist")).readlines()]
        raw_vcfs = [x for x in raw_vcfs if re.sub("_v[1-9]", "", x.split("/")[-1]) in blacklist] #TODO: remove the sub once the bam names are fixed
    if config.get("whitelist"):
        whitelist = [x.strip() for x in open(config.get("whitelist")).readlines()]
        raw_vcfs = [x for x in raw_vcfs if re.sub("_v[1-9]", "", x.split("/")[-1]) in whitelist] #TODO: remove the sub once the bam names are fixed

logging.debug("%s input files left after blacklist/whitelist filtering." % len(raw_vcfs))

all_maf_name = config["raw_vcf_folder"] + "all_samples_%s_S%s.maf" % (config["name"], len(raw_vcfs))
all_annotated_name = config["raw_vcf_folder"] + "all_samples_%s_S%s.tsv" % (config["name"], len(raw_vcfs))

@follows(clean_folders_of_intermediate_files)
@transform(raw_vcfs, suffix(".vcf"), "_annotated.tsv")
def annotate_raw_vcfs(infile, outfile):
    sample = SA.SampleAnnotation(infile, target_folder= "", run_tools=True)
    sample.print_rows()

@follows(annotate_raw_vcfs)
@merge(annotate_raw_vcfs, all_annotated_name)
def unite_annotated_raw_vcfs(samples, output):
    SA.annotate_all_samples_as_one(filenames=samples, outfile=all_annotated_name)

###########################################################
# VCF2MAF for unfiltered vcfs

@follows(unite_annotated_raw_vcfs)
@transform(raw_vcfs, suffix(".vcf"), ".maf")
def transform_raw_vcfs_to_mafs(infile, outfile):
    vcf2maf.run_vcf2maf(infile)

@follows(transform_raw_vcfs_to_mafs)
@merge(transform_raw_vcfs_to_mafs, all_maf_name)
def unite_mafs(samples, output):
    logging.debug("Uniting:", samples)
    vcf2maf.unite_mafs(samples, output)

########################################################################################
# Filter the raw vcfs ##################################################################
########################################################################################

@follows(unite_annotated_raw_vcfs)
@follows(unite_mafs, mkdir(config["vcf_folder"]))
@transform(raw_vcfs, regex(r'.+\/(.+)'), r"%s\1" % config["vcf_folder"], config["vcf_folder"])
def filter_vcf(infile, outfile, outdir):
    
    logging.debug("Filtering %s to %s" % (infile, outfile))
    output = open(outfile, "w")
    
    min_cov = int(config.get("min_cov", 0))
    min_varfreq = float(config.get("min_varfreq", 0))
    min_qual = int(config.get("min_qual", 0))
    
    try:
        vcf_type = config["vcf_type"]
    except:
        print "Sorry, you need to specify a vcf_type in the config file, e.g. 'illumina_strelka' or 'iontorrent'."
            
    with open(infile) as handle:
        variants = 0
        
        for row in handle:
            if row[0] == "#":
                output.write(row)
            else:
                
                if not config.get("no_filtering"):
                    
                    if vcf_type == "iontorrent":
                        if "CNV" not in row and "Non-Confident" not in row: # low confidence due to low Normal coverage
                            var = AP.VCFrow(row)
                            fao, fdp = float(var.values.get("FAO", 0)), float(var.values.get("FDP", 0))
                            
                            if min_cov < fao+fdp  and min_varfreq < fao/fdp:
                                variants += 1
                                output.write(row)
                            
                    elif vcf_type == "illumina_strelka":
                        var = AP.VCFrow(row)
                        
                        dp = int(var.values.get("DP", 0))
                        qss = int(var.values.get("QSS", 0))
                        
                        if min_cov < dp and min_qual < qss:
                            variants += 1
                            output.write(row)
                            
                    else:
                        logging.critical( "Sorry, the vcf_type specified in the config file must be known by the main pipeline script." )
                        logging.critical( "Please define it in the filter_vcf() function of the pipeline script." )
                        logging.critical( "Unknown vcf_type: %s" % vcf_type )
                        
                else:
                    variants += 1
                    output.write(row)
    
    output.close()
    #
    #if variants == 0:
    #    os.remove(outdir + outfile)

####################################################################
# Fetch annotation for filtered vcfs from unfiltered vcfs.

all_filtered_maf_name = config["vcf_folder"] + "all_samples_%s_S%s.maf" % (config["name"], len(raw_vcfs))
all_filtered_annotated_name = config["vcf_folder"] + "all_samples_%s_S%s.tsv" % (config["name"], len(raw_vcfs))

@follows(filter_vcf)
@transform(filter_vcf, suffix(".vcf"), "_annotated.tsv")
def fetch_annotation_for_filtered_vcfs(infile, outfile):
    sample = infile[:-4].split("/")[-1]
    raw_annotated = config["raw_vcf_folder"] + sample + "_annotated.tsv"
    
    accepted_chromposes = []
    with open(infile) as filtered: 
        for row in filtered:
            var = AP.VCFrow(row)
            if var.chrompos:
                accepted_chromposes.append(var.chrompos)
    accepted_chromposes = set(accepted_chromposes)
    
    out = open(outfile, "w")
    accepted = 0
    with open(raw_annotated) as handle:
        for row in handle:
            if row[0] == "#" or "chrom" in row:
                out.write(row)
            else:
                f = row.split("\t")
                chrompos = f[1] +"\t"+ f[2]
                if chrompos in accepted_chromposes:
                    accepted += 1 
                    out.write(row)
    out.close()
    logging.debug("Fetched %s of %s annotated rows for %s from %s." % (accepted, len(accepted_chromposes), infile, outfile))


@follows(fetch_annotation_for_filtered_vcfs)
@merge(fetch_annotation_for_filtered_vcfs, all_filtered_annotated_name)
def unite_annotated_filtered_vcfs(samples, output):
    logging.debug("Uniting %s annotated samples starting from %s into %s." % (len(samples), samples[0], output))
    SA.annotate_all_samples_as_one(filenames=samples, outfile=all_filtered_annotated_name)

####################################################################
# VCF2MAF for filtered vcfs

@follows(unite_annotated_filtered_vcfs)
@transform(filter_vcf, suffix(".vcf"), ".maf")
def transform_filtered_vcfs_to_mafs(infile, outfile):
    logging.debug("Starting vcf2maf from %s to %s." % (infile, outfile))
    vcf2maf.run_vcf2maf(infile)

@follows(transform_filtered_vcfs_to_mafs)
@merge(transform_filtered_vcfs_to_mafs, all_filtered_maf_name)
def unite_filtered_mafs(samples, output):
    vcf2maf.unite_mafs(samples, output)


########################################################################################
# Functional analysis programs #########################################################
########################################################################################

analysis_folder = config["analysis_folder"]
music_run_flag = config["analysis_folder"] + "music_run_flag" # so ruffus has a way of knowing when Music has run
mutsig_run_flag = config["analysis_folder"] + "mutsig_run_flag" # so ruffus has a way of knowing when MutSigCV has run
mutsig_rename_flag = config["analysis_folder"] + "mutsig_rename_flag" # so ruffus has a way of knowing the files were renamed
mutsig_output = config["analysis_folder"] + "/output/mutsigcv"

@active_if(config.get("bam_folder"))
@follows(unite_filtered_mafs)
@follows(unite_annotated_filtered_vcfs)
@merge([unite_filtered_mafs, unite_annotated_filtered_vcfs], music_run_flag)
def run_music(input_files, music_run_flag):
    logging.debug("Starting MuSiC from %s." % input_files)
    music.prepare_input_and_run_music(input_maf=input_files[0], annotationfile=input_files[1], config=config, runflag=music_run_flag)

@follows(unite_filtered_mafs)
@transform(unite_filtered_mafs, regex(r'(.+\/).+'), mutsig_run_flag)
def run_mutsig(input_maf, mutsig_run_flag):
    logging.debug("Starting MutSigCV from %s." % input_maf)
    mutsig.run_mutsig(input_maf, mutsig_output)
    open(mutsig_run_flag, "w")

@follows(run_mutsig)
@transform(run_mutsig, regex(r'(.+\/).+'), mutsig_rename_flag)
def rename_mutsig_out(infile, mutsig_rename_flag):
    output_folder = config["analysis_folder"] + "/output/"
    for filename in os.listdir(output_folder):
        
        pipeline_identifier = project_name + "_S " + str(len(raw_vcfs))
        
        if os.path.isdir(output_folder+filename) or filename.startswith(pipeline_identifier): continue        
        new_name = pipeline_identifier + "_" + filename
        os.rename(output_folder+filename, output_folder+new_name)
    open(mutsig_rename_flag, "w")

        
########################################################################################
########################################################################################
########################################################################################

last_steps = [run_music, rename_mutsig_out]

pipeline_printout(sys.stdout, last_steps)
pipeline_printout_graph ( open("pipeline.png", "w"),"png", last_steps, no_key_legend=True)
pipeline_run(multiprocess = config.get("cpus", 1))