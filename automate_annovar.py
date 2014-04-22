# Get a filename as an input
# Convert the vcf to an AnnoVar input file
# Run AnnoVar
# Parse and return the output

import sys
import os 
import subprocess

def run_annovar(target_folder, target_vcf):
    outfile_name = target_vcf[:-4]
    convert_cmd = """perl /home/andreas/bioinfo/projects/wtc_quasar_analysis/scripts/annovar/convert2annovar.pl --format vcf4 %s > %s.annoin""" % (target_folder+target_vcf, target_folder+outfile_name)
    
    if "vcf4old" in convert_cmd:
        print "Careful! Annoin generation uses --format vcf4old flag!"
     
    subprocess.call(convert_cmd, shell=True, stdout = open("log_out.txt", "wa"), stderr = open("log_err.txt", "wa"))
    #print convert_cmd 
    
    annotate_cmd = """perl /home/andreas/bioinfo/projects/wtc_quasar_analysis/scripts/annovar/table_annovar.pl %s.annoin /home/andreas/bioinfo/projects/wtc_quasar_analysis/scripts/annovar/humandb/ -buildver hg19 -out %s -remove -protocol refGene,phastConsElements46way,genomicSuperDups,esp6500si_all,1000g2012apr_all,1000g2012apr_eur,1000g2012apr_afr,1000g2012apr_amr,1000g2012apr_asn,snp135,ljb2_all,cosmic67 -operation g,r,r,f,f,f,f,f,f,f,f,f -nastring - """ % (target_folder+outfile_name, target_folder+outfile_name)
    
    #print annotate_cmd 
    subprocess.call(annotate_cmd, shell=True, stdout = open("log_out.txt", "wa"), stderr = open("log_err.txt", "wa"))
    
    annovar_output = target_folder + outfile_name + ".hg19_multianno.txt"
    
    return annovar_output 
    
##############################################################################################

if __name__ == '__main__':

    target_folder = os.getcwd() + "/"
    target_vcf = sys.argv[1]
    
    run_annovar(target_folder, target_vcf)


# annovar_output = target_folder + outfile_name + ".hg19_multianno.txt"
# 
# rows = {}
# 
# with open(annovar_output) as output:
#     for row in output:
# #         print row
#         
#         if "Start" in row:
#             header = row.strip("\n").split("\t")
#             
#         else:
#             f = row.strip("\n").split("\t")
#             chrompos = f[0] +"\t"+ f[1]
#             
#             rows[chrompos] = f
#             
#             for field in f:
#                 if "COSM" in field:
#                     print field
#                     cosmic_no = field.split(";")[0]
#                     
#                     entries =  field.split(";")[1:]
#                     for entry in entries:
#                         occurence = entry.split("(")[0]
#                         tissue = entry.split("(")[1][:-2]
# 
#                         print cosmic_no, occurence, tissue