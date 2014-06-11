# Get a filename as an input
# Run snpeff

import sys
import os 
import subprocess
from ToolConfig import snpeff_jar, snpeff_config_file
import logging

def run_snpeff(target_folder, target_vcf):
    
    outfile_name = target_vcf[:-4]
    
    snpeff_input = target_folder + target_vcf
    snpeff_output = target_folder + outfile_name + ".snpeff.vcf"
  
    annotate_cmd = """java -Xmx2g -jar %s hg19 -t -hgvs -c %s -v %s > %s""" % (snpeff_jar, snpeff_config_file, snpeff_input, snpeff_output)
    
    logging.debug( annotate_cmd ) 
    subprocess.call(annotate_cmd, shell=True, stdout = open("log_out.txt", "wa"), stderr = open("log_err.txt", "wa"))
    
    return snpeff_output 
    
##############################################################################################

if __name__ == '__main__':

    target_folder = os.getcwd() + "/"
    target_vcf = sys.argv[1]
    
    run_snpeff(target_folder, target_vcf)
    
    
