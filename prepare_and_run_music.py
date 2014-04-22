# parse the available bams, mafs and coverage files
# find the minimal list of samples that have an entry in both and can thus be used for MUsic
# create missing covg files if necessary
# trim down all list to this set
#
# run Genome MuSiC
#

import sys
import os
import pandas as pd
import subprocess
import multiprocessing as mp
import re
import pprint

def prepare_music_input(config=None, input_maf=None, annotationfile=None, runflag=None):
    
    pprint.pprint(config)
    
    bam_folder = config["bam_folder"]
    
    analysis_folder = config["analysis_folder"]
    coverage_folder =  analysis_folder + "output/roi_covgs"
    output_folder = analysis_folder + "output/"
    input_folder = analysis_folder + "input/"
    bedfile = config["bed"]
    reffile = config["ref"]
    
    kegg_pathwayfile = input_folder + "KEGG_052412.txt"
    react_pathwayfile = input_folder + "reactome_pathway_file"
    
    #root = "/home/andreas/bioinfo/projects/wtc_quasar_analysis/data/"
    #bam_folder = root + "quasar_bams/"
    #input_maf = root + "quasar_vcfs_default_parameters_QC2/all_samples_QC2.maf"
    #coverage_folder = root + "music_QC2/output/roi_covgs/"
    #output_folder = root + "music_QC2/output/"
    #input_folder = root + "music_QC2/input/"
    #bedfile = root + "music_QC2/input/roi_148gene_panel_HP.bed"
    #reffile = "/home/andreas/bioinfo/core/general/data/hg19.fasta"
    #annotationfile = root + "quasar_vcfs_default_parameters_QC2/all_samples_default_parameters_QC2.tsv"
    #kegg_pathwayfile = root + "music_QC2/input/pathways/KEGG_052412.txt"
    #react_pathwayfile = root + "music_QC2/input/pathways/reactome_pathway_file"
    
    cpus = 1
   
    ######################################################################################
    # fetch bams/maf and parse bams into dict of paths
    
    bam_files = [x for x in os.listdir(bam_folder)]
    bam_samples = set([unify_sample_name(x) for x in bam_files if x.endswith("T.bam")])
    bam_paths = {bam:dict(Normal = None, Tumor = None) for bam in bam_samples}
       
    for bam in bam_files:
        if not bam_paths.get(unify_sample_name(bam)):
            continue
        
        if bam.endswith("N.bam"):
            bam_paths[unify_sample_name(bam)]["Normal"] = bam_folder+bam
        elif bam.endswith("T.bam"):
            bam_paths[unify_sample_name(bam)]["Tumor"] = bam_folder+bam
    
    maf_df = pd.read_csv(input_maf, sep = "\t")
    maf_df = maf_df.dropna(subset=["Tumor_Sample_Barcode"])
    maf_samples = set([unify_sample_name(x) for x in maf_df.Tumor_Sample_Barcode.dropna().unique()])

    ######################################################################################
    # process
    
    music_input_samples = bam_samples & maf_samples
    music_input_samples = set([x for x in music_input_samples if bam_paths.get(x)]) # reduce to bams with dict entry
    print "Number of samples present in both bamfolder and maf-file:", len(music_input_samples)
    
    bams_output_file = input_folder + "music_bam_path_list_%s.txt" % len(music_input_samples)
    
    coverage_samples = set([unify_sample_name(x) for x in os.listdir(coverage_folder)])
    coverage_2do_samples = music_input_samples - coverage_samples 
    
    if coverage_2do_samples:
        print "Found %s bams that need coverage calculation, starting..." % len(coverage_2do_samples)
        run_multiple_calc_covg(ref=reffile, bed = bedfile, bams = coverage_2do_samples, bam_dict = bam_paths, outfolder = coverage_folder, cpus=cpus)
        
    print "Uniting covg files..."
    bams_path_file = create_music_bam_path_file(music_input_samples, bam_paths, bams_output_file)
    unite_calc_covg(ref=reffile, bed=bedfile, bams=bams_path_file, outfolder=output_folder)

    ######################################################################################
    # delete covg files for which no maf entry is found
    # as e.g. some covg file copied over from another project might have been excluded 
    
    # seems unnecessary, no error if skipped
    
    #for covg_file in os.listdir(coverage_folder):
    #    if unify_sample_name(covg_file) not in maf_samples:
    #        
    #        print "Deleting", coverage_folder+covg_file
    #        os.remove(coverage_folder+covg_file)
    
    ######################################################################################
    # create output
    
    print "Printing bam list for MuSiC"
    bams_path_file = create_music_bam_path_file(music_input_samples, bam_paths, bams_output_file)
    
    print "Printing maf input for MuSiC"
        
    modified_maf = input_folder + "music_input_%s.maf" % len(music_input_samples)
    create_modified_maf(bams=music_input_samples, maf=input_maf, outfile=modified_maf, annotationfile=annotationfile)
    #maf_df["Tumor_Sample_Barcode"] = maf_df.apply(lambda x: unify_sample_name(x["Tumor_Sample_Barcode"]), axis=1)
    #maf_df["Matched_Norm_Sample_Barcode"] = maf_df.apply(lambda x: unify_sample_name(x["Matched_Norm_Sample_Barcode"])+"_N", axis=1)
    #mask = maf_df.apply(lambda x: unify_sample_name(x["Tumor_Sample_Barcode"]) in music_input_samples, axis=1)
    #maf_df[mask].to_csv(modified_maf, sep="\t")
    
    #######################################################################################
    ## MuSiC
    
    print "#" * 50
    print "################ Starting MuSiC ##################"
    print "#" * 50
    
    run_bmr(ref=reffile, bed=bedfile, bams=bams_path_file, outfolder=output_folder, maf=modified_maf)
    run_smg(outfolder=output_folder)
    run_mutation_relation(outfolder=output_folder, maf=modified_maf, bams=bams_path_file)
    
    for pathway in [kegg_pathwayfile, react_pathwayfile]:
        run_pathways(outfolder=output_folder, maf=modified_maf, bams=bams_path_file, pathways=pathway)
    
    run_proximity(outfolder=output_folder, maf=modified_maf)
    run_omim(outfolder=output_folder, maf=modified_maf)
    
    #######################################################################################
    # create success file for ruffus
    
    flag = open(runflag, "w")
    flag.close()
    print "MuSiC run succesfully!"
    
######################################################################################

def create_modified_maf(bams=None, maf=None, outfile=None, annotationfile=None):
    
    ##############################################################
    # parse annotation into dict to add the 3 required Music columns later
    
    anno = {}
    with open(annotationfile) as handle:
        for row in handle:
            if "sample" in row: continue
            #print row
            f = row.split()
            chrompos_refalt = "\t".join([f[1], f[2], f[4], f[5]])
            info = f[53]
        
            if info != "-":
                #print info
                info = info.split(",")[0].split(":")
                
                if len(info) > 4:
                    anno[chrompos_refalt] = (info[1], info[3],info[4]) 
        
    ##############################################################
    
    output = open(outfile, "w")
        
    with open(maf) as handle:
        for row in handle:
            
            if "Tumor_Sample_Barcode" in row:
                row = row.strip() + "\ttranscript_name\tamino_acid_change\tc_position"
                output.write(row+"\n")
                continue
            
            elif row[0] == "#":
                #output.write(row)
                continue

            f = row.split()
            tumor = f[13]
            chrompos_refalt = "\t".join([f[4], f[5], f[10], f[11]])
            
            if unify_sample_name(tumor) in bams:
                extra_info = anno.get(chrompos_refalt, ("-", "-", "-"))
                row = row.strip() +"\t"+ "\t".join(extra_info)
                row = re.sub("_v[1-9]", "", row)
                output.write(row + "\n")
    
    output.close()

def unify_sample_name(raw):
    """
    Modifiy the sample names to be comparable
    """
    sample = raw.split(".")[0]
    sample = "_".join(sample.split("_")[:2])
    
    return sample

def create_music_bam_path_file(bam_list, bam_dict, output_file):
    
    with open(output_file, "w") as output:
        for bam in bam_list:
            result = "\t".join([bam, bam_dict[bam]["Normal"], bam_dict[bam]["Tumor"]])
            output.write(result + "\n")
    
    return output_file

def unite_calc_covg(ref=None, bed = None, bams = None, outfolder = None):
    unite_cmd = "genome music bmr calc-covg --roi-file %s --reference-sequence %s --bam-list %s --output-dir %s" % (bed, ref, bams, outfolder)
    #print unite_cmd
    out = open("music_log.txt", "wa")
    subprocess.call(unite_cmd, shell=True, stdout=out, stderr=out)
    out.close()
    
def run_multiple_calc_covg(bed = None, ref=None, bams = None, bam_dict = None, outfolder = None, cpus=None):
    
    if cpus > 1:
        pool = mp.Pool(cpus)
    
    for bam in bams:
        nbam = bam_dict[bam]["Normal"]
        tbam = bam_dict[bam]["Tumor"]
        out = outfolder + bam + ".covg"
    
        calc_cmd = "calcRoiCovg %s %s %s %s %s 6 8 20" % (nbam, tbam, bed, ref, out)

        if cpus == 1:
            execute(calc_cmd, bam)
        else:
            pool.apply_async(execute, args = [calc_cmd, bam])
    
    if cpus > 1:
        pool.close()
        pool.join()
        
def execute(cmd, sample=None):
    if sample:
        print "Starting process for", sample
        
    subprocess.call(cmd, shell=True)
        
def run_bmr(bams=None, maf=None, outfolder=None, ref=None, bed=None):
    print "#" * 50
    print "Starting background mutation rate calculation (bmr)."
    
    bmr_cmd = """genome music bmr calc-bmr \
    --bam-list %s \
    --maf-file %s \
    --output-dir %s \
    --reference-sequence %s \
    --roi-file %s \
    --genes-to-ignore TP53,APC """ % (bams, maf, outfolder, ref, bed)
    
    out = open("music_log.txt", "wa")
    subprocess.call(bmr_cmd, shell=True, stdout=out, stderr=out)
    out.close()
    
def run_smg(outfolder=None):
    print "#" * 50
    print "Starting calculation of sign. mutated genes (smg)."
    
    gene_mrs = outfolder + "gene_mrs"
    out = outfolder + "smgs"

    smg_cmd = """genome music smg \
          --gene-mr-file %s \
          --output-file %s """ % (gene_mrs, out)
    
    #print smg_cmd
    out = open("music_log.txt", "wa")
    subprocess.call(smg_cmd, shell=True, stdout=out, stderr=out)
    out.close()
    
def run_mutation_relation(outfolder=None, maf=None, bams=None):
    print "#" * 50
    print "Starting permutation of related mutations."
    
    outfile = outfolder + "mutation_relation.csv"
    matrix = outfolder + "mutation_relation_matrix.csv"
    
    relation_cmd = """genome music mutation-relation \
        --bam-list %s \
        --maf-file %s \
        --permutations 1000 \
        --mutation-matrix-file %s\
        --output-file %s""" % (bams, maf, matrix, outfile) 
    
    #print relation_cmd
    out = open("music_log.txt", "wa")
    subprocess.call(relation_cmd, shell=True, stdout=out, stderr=out)
    out.close()

def run_pathways(outfolder=None, maf=None, bams=None, pathways=None):
    
    pathway_name = pathways.split("/")[-1].split(".")[0]
    gene_covg_dir = outfolder + "gene_covgs"
    out = outfolder + "sm_pathways_" + pathway_name 

    print "#" * 50
    print "Starting pathway analysis for %s." % (pathway_name)

    path_cmd = """genome music path-scan \
        --bam-list %s \
        --gene-covg-dir %s \
        --maf-file %s \
        --output-file %s \
        --pathway-file %s \
        --genes-to-ignore TP53 \
        --bmr 8.7E-07""" % (bams, gene_covg_dir, maf, out, pathways)
        
    out = open("music_log.txt", "wa")
    subprocess.call(path_cmd, shell=True, stdout=out, stderr=out)
    out.close() 

def run_proximity(outfolder=None, maf=None):
    
    print "#" * 50
    print "Starting proximity analysis." 

    prox_cmd = """genome music proximity \
        --maf-file %s \
        --output-dir %s \
        --max-proximity 150""" % (maf, outfolder)
    
    out = open("music_log.txt", "wa")
    subprocess.call(prox_cmd, shell=True, stdout=out, stderr=out)
    out.close()

def run_omim(outfolder=None, maf=None):
    
    print "#" * 50
    print "Starting OMIM annotation." 

    outfile = outfolder + "omim_annotation.tsv"

    omim_cmd = """genome music cosmic-omim \
        --maf-file %s \
        --output-file %s \
        --no-verbose""" % (maf, outfile)
    
    #print omim_cmd    
    out = open("music_log.txt", "wa")
    subprocess.call(omim_cmd, shell=True, stdout=out, stderr=out)
    out.close()

##############################################################################################

if __name__ == '__main__':
    main()