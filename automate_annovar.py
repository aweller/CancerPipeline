# Get a filename as an input
# Convert the vcf to an AnnoVar input file
# Run AnnoVar
# Parse and return the output

import sys
import os 
import subprocess
from ToolConfig import convert2annovar_script, table_annovar_script, annovar_db_folder
import logging

def run_annovar(target_folder, target_vcf):
    outfile_name = target_vcf[:-4]
    convert_cmd = """perl %s --format vcf4 %s > %s.annoin""" % (convert2annovar_script, target_folder+target_vcf, target_folder+outfile_name)
    
    if "vcf4old" in convert_cmd:
        logging.warning( "Careful! Annoin generation uses --format vcf4old flag!" )
     
    subprocess.call(convert_cmd, shell=True, stdout = open("log_out.txt", "wa"), stderr = open("log_err.txt", "wa"))
    logging.debug( convert_cmd )
    
    #annotate_cmd = """perl %s %s.annoin %s -buildver hg19 -out %s -remove -protocol %s -operation %s -nastring - """ % (table_annovar_script, target_folder+outfile_name, annovar_db_folder, target_folder+outfile_name, annovar_protocols, annovar_operation)

    annotate_cmd = """perl %s %s.annoin %s -buildver hg19 -out %s -remove \
                   -protocol refGene,phastConsElements46way,genomicSuperDups,esp6500si_all,1000g2012apr_all,1000g2012apr_eur,1000g2012apr_afr,1000g2012apr_amr,1000g2012apr_asn,snp135,ljb2_all,cosmic67\
                   -operation g,r,r,f,f,f,f,f,f,f,f,f\
                   -nastring - """ % (table_annovar_script, target_folder+outfile_name, annovar_db_folder, target_folder+outfile_name)

    logging.debug( annotate_cmd )
    result = subprocess.call(annotate_cmd, shell=True, stdout = open("log_out.txt", "wa"), stderr = open("log_err.txt", "wa"))
    
    if result > 0:
        logging.critical( annotate_cmd )
        logging.critical( "Annovar failed. Is this a working command?" )
        sys.exit()
        
    annovar_output = target_folder + outfile_name + ".hg19_multianno.txt"
    
    return annovar_output 
    
##############################################################################################

if __name__ == '__main__':
    
    logger=logging.getLogger()
    logger.setLevel(logging.DEBUG)

    target_folder = os.getcwd() + "/"
    target_vcf = sys.argv[1]
    
    run_annovar(target_folder, target_vcf)