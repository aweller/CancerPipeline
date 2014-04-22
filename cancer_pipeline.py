from ruffus import *
import sys
import os
import AnnotationRowParsers as AP
import SampleAnnotation as SA
import automate_vcf2maf as vcf2maf
from cancer_pipeline_scripts import *
import prepare_and_run_music as music
import automate_mutsig as mutsig
import pipeline_config as config
import pprint
import shutil
import logging

# TODO
# - implement a fast elegant way that filters the vcfs and then creates maf from the old mafs.
# - SQL integration
# - proper debugging output
# - use sortedDict in the annotated output

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

########################################################################################
# Prepare the unfiltered vcfs ##########################################################
########################################################################################

###########################################################
# Annotation

raw_vcfs = [config["raw_vcf_folder"]+x for x in os.listdir(config["raw_vcf_folder"])
                  if x.endswith(".vcf") and
                  "QC" not in x and
                  "all" not in x and
                  "anno" not in x]

if config.get("blacklist"):
    raw_vcfs = [x for x in raw_vcfs if x not in config.get("blacklist")]

logging.debug("Found %s raw input files:" % len(raw_vcfs))
logging.debug(raw_vcfs)

#all_maf_name = config["raw_vcf_folder"] + "all_samples_%s_S%s.maf" % (config["name"], len(raw_vcfs))
#all_annotated_name = config["raw_vcf_folder"] + "all_samples_%s_S%s.tsv" % (config["name"], len(raw_vcfs))
#
#@transform(raw_vcfs, suffix(".vcf"), "_annotated.tsv")
#def annotate_raw_vcfs(infile, outfile):
#    sample = SA.SampleAnnotation(infile, target_folder= "", run_tools=True)
#    sample.print_rows()
#
#@follows(annotate_raw_vcfs)
#@merge(annotate_raw_vcfs, all_annotated_name)
#def unite_annotated_raw_vcfs(samples, output):
#    SA.annotate_all_samples_as_one(filenames=samples, outfile=all_annotated_name)
#
############################################################
## VCF2MAF for unfiltered vcfs
#
#@transform(raw_vcfs, suffix(".vcf"), ".maf")
#def transform_raw_vcfs_to_mafs(infile, outfile):
#    vcf2maf.run_vcf2maf(infile)
#
#@follows(transform_raw_vcfs_to_mafs)
#@merge(transform_raw_vcfs_to_mafs, all_maf_name)
#def unite_mafs(samples, output):
#    vcf2maf.unite_mafs(samples, output)

########################################################################################
# Filter the raw vcfs ##################################################################
########################################################################################

#@follows(unite_annotated_raw_vcfs)
#@follows(unite_mafs, mkdir(config["vcf_folder"]))

@transform(raw_vcfs, regex(r'.+\/(.+)'), r"%s\1" % config["vcf_folder"], config["vcf_folder"], int(config["max_snv"]), int(config["min_cov"]), float(config["min_varfreq"]))
def filter_vcf(infile, outfile, outdir, max_snv, min_cov, min_varfreq):
    
    logging.debug("Filtering %s to %s" % (infile, outfile))

    #outfile = outfile.split("/")[-1]    
    deaminated = False
    #snv_no = len([x for x in open(infile).readlines() if "snp" in x])
    #if snv_no > max_snv:
    #    deaminated = True

    output = open(outfile, "w")
            
    if not deaminated:
        with open(infile) as handle:
            variants = 0
            
            for row in handle:
                if row[0] == "#":
                    output.write(row)
                else:
                    if "CNV" not in row and "Coverage" not in row: # low confidence due to low Normal coverage
                        var = AP.VCFrow(row)
                        fao, fdp = float(var.values["FAO"]), float(var.values["FDP"])
                        if min_cov < fao+fdp  and min_varfreq < fao/fdp:
                            variants += 1
                            output.write(row)                     
    
    output.close()
    #
    #if variants == 0:
    #    os.remove(outdir + outfile)

####################################################################
# repeat the maf creation for the filtered vcf. This is a slow and silly way, but quicker right now.

all_filtered_maf_name = config["vcf_folder"] + "all_samples_%s_S%s.maf" % (config["name"], len(raw_vcfs))
all_filtered_annotated_name = config["vcf_folder"] + "all_samples_%s_S%s.tsv" % (config["name"], len(raw_vcfs))

@follows(filter_vcf)
@transform(filter_vcf, suffix(".vcf"), "_annotated.tsv")
def annotate_filtered_vcfs(infile, outfile):
    sample = SA.SampleAnnotation(infile, target_folder= "", run_tools=True)
    sample.print_rows()

@follows(annotate_filtered_vcfs)
@merge(filter_vcf, all_filtered_annotated_name)
def unite_annotated_filtered_vcfs(samples, output):
    SA.annotate_all_samples_as_one(filenames=samples, outfile=all_filtered_annotated_name)

####################################################################
# VCF2MAF for filtered vcfs

@follows(filter_vcf)
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
mutsig_output = config["analysis_folder"] + "/output/mutsigcv_"

@active_if(config.get("bam_folder"))
@follows(unite_filtered_mafs)
@follows(unite_annotated_filtered_vcfs,
         mkdir([analysis_folder,
                analysis_folder+"/input",
                analysis_folder+"/output",
                analysis_folder+"/output/roi_covgs",
                analysis_folder+"/output/gene_covgs"]))
@merge([unite_filtered_mafs, unite_annotated_filtered_vcfs], music_run_flag)
def run_music(input_files, music_run_flag):
    logging.debug("Starting MuSiC from %s." % input_files)
    music.prepare_music_input(input_maf=input_files[0], annotationfile=input_files[1], config=config, runflag=music_run_flag)

@follows(unite_filtered_mafs)
@transform(unite_filtered_mafs, regex(r'(.+\/).+'), r"\1%s" % "mutsig_run_flag")
def run_mutsig(input_maf, mutsig_run_flag):
    logging.debug("Starting MutSigCV from %s." % input_maf)
    mutsig.run_mutsig(input_maf, mutsig_output)

########################################################################################
########################################################################################
########################################################################################

last_steps = [run_music, run_mutsig]

pipeline_printout(sys.stdout, last_steps)
pipeline_printout_graph ( open("pipeline.png", "w"),"png", last_steps, no_key_legend=True)
pipeline_run(multiprocess = 4)