import os
import subprocess
import sys

def run_mutsig(input_maf, output_file):
    
    db = "/home/andreas/bioinfo/core/mutsig/database"
    output_prefix = output_file    
    matlab_dir = "/home/andreas/bin/matlab2/MATLAB_Compiler_Runtime/v81/"
    
    mutseq_cmd = """/home/andreas/bioinfo/core/mutsig/MutSigCV_1.4/run_MutSigCV.sh %s %s %s/exome_full192.coverage.txt %s/gene.covariates.txt %s %s/mutation_type_dictionary_file.txt %s/chr_files_hg19""" % (matlab_dir, input_maf, db, db, output_prefix, db, db)

    print mutseq_cmd
    subprocess.call(mutseq_cmd, shell = True)
    
###########################################################

if __name__ == '__main__':
   input_maf = sys.argv[1]
   output_file = sys.argv[2]
   run_mutsig(input_maf, output_file)