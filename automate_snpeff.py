# Get a filename as an input
# Run snpEff

import sys
import os 
import subprocess

def run_snpeff(target_folder, target_vcf):
    
    outfile_name = target_vcf[:-4]
    
    snpEff_path = "/home/andreas/bioinfo/projects/wtc_quasar_analysis/scripts/snpeff"
    snpEff_input = target_folder + target_vcf
    snpEff_output = target_folder + outfile_name + ".snpeff.vcf"
  
    annotate_cmd = """java -Xmx2g -jar %s/snpEff.jar hg19 -t -hgvs -c %s/snpEff.config -v %s > %s""" % (snpEff_path, snpEff_path, snpEff_input, snpEff_output)
    
    #print annotate_cmd 
    subprocess.call(annotate_cmd, shell=True, stdout = open("log_out.txt", "wa"), stderr = open("log_err.txt", "wa"))
    
    return snpEff_output 
    
##############################################################################################

if __name__ == '__main__':

    target_folder = os.getcwd() + "/"
    target_vcf = sys.argv[1]
    
    run_snpEff(target_folder, target_vcf, outfile_name)
    
    
