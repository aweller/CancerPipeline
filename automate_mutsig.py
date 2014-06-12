import os
import subprocess
import sys
from ToolConfig import mutsig_db, mutsig_matlab_folder, mutsig_run_script
import logging

def run_mutsig(input_maf, output_file):
    
    output_prefix = output_file    
    
    mutseq_cmd = """%s %s %s %s/exome_full192.coverage.txt %s/gene.covariates.txt \
                 %s %s/mutation_type_dictionary_file.txt %s/chr_files_hg19""" % (mutsig_run_script, mutsig_matlab_folder, input_maf,
                                                                                 mutsig_db, mutsig_db, output_prefix, mutsig_db, mutsig_db)
    logging.debug( mutseq_cmd )
    result = subprocess.call(mutseq_cmd, shell = True)
    
    if result > 0:
        logging.critical( mutseq_cmd )
        logging.critical( "MutSigCV failed. Is this a working command?" )
        sys.exit()
    
###########################################################

if __name__ == '__main__':
   input_maf = sys.argv[1]
   output_file = sys.argv[2]
   run_mutsig(input_maf, output_file)