from ruffus import *
import sys
import os

starting_files = [x for x in os.listdir("./")
                  if x.endswith(".vcf") and
                  "QC" not in x]

@transform(starting_files, suffix(".vcf"), "_QC.vcf")
def filter_vcf(infile, outfile):
    #infile = kwargs["input"]
    #outfile = kwargs["output"]
    
    print "Filtering vcf."
    
    out = open(outfile, "w")
    out.close()
    
    return outfile

@transform(filter_vcf, suffix("_QC.vcf"), "_QC.maf")
def vcf2maf(infile, outfile):
    #infile = kwargs["input"]
    #outfile = kwargs["output"]
    
    print "Converting vcf to maf."
    
    out = open(outfile, "w")
    out.close()
    
    return outfile


@merge(vcf2maf, 'all.summary')
def unite_all(targets, outfile):
    
    print "Uniting %s targets." % len(targets)
    
    out = open(outfile, "w")
    out.close()

pipeline_printout(sys.stdout, [unite_all], verbose = 3)
pipeline_run()
#pipeline_printout_graph ( open("simple_tutorial_stage5_before.png", "w"),"png",[unite_all],minimal_key_legend=True)