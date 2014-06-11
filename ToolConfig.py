#############################################################################################################
# vcf2maf
# This is a script from Github which you can fetch with
# git clone https://github.com/ckandoth/vcf2maf.git

vcf2maf_script = "/home/andreas/bioinfo/core/vcf2maf_perl/vcf2maf-master/vcf2maf.pl"

#############################################################################################################
# SNPeff

snpeff_jar = "/home/andreas/bioinfo/projects/wtc_quasar_analysis/scripts/snpeff/snpEff.jar"
snpeff_config_file = "/home/andreas/bioinfo/projects/wtc_quasar_analysis/scripts/snpeff/snpEff.config"

#############################################################################################################
# Annovar

convert2annovar_script = "/home/andreas/bioinfo/projects/wtc_quasar_analysis/scripts/annovar/convert2annovar.pl"
table_annovar_script = "/home/andreas/bioinfo/projects/wtc_quasar_analysis/scripts/annovar/table_annovar.pl"
annovar_db_folder = "/home/andreas/bioinfo/projects/wtc_quasar_analysis/scripts/annovar/humandb/"

# the protocols specify from which information sources Annovar should fetch (see Annovar documentation)
annovar_protocols = """refGene,phastConsElements46way,genomicSuperDups,esp6500si_all,
    1000g2012apr_all,1000g2012apr_eur, 1000g2012apr_afr,1000g2012apr_amr,1000g2012apr_asn,snp135,ljb2_all,cosmic67"""
# the operation must contain one field for each field in the protocols and specifies the data type (see Annovar documentation)
annovar_operation = "g,r,r,f,f,f,f,f,f,f,f,f"

#############################################################################################################
# MutSigCV
#
# The matlab support needed for MutSigCV can be a pain to install, if it doesn't work try a different version (e.g. v81 instead of v83)

mutsig_db = "/home/andreas/bioinfo/core/mutsig/database"
mutsig_matlab_folder = "/home/andreas/bin/matlab2/MATLAB_Compiler_Runtime/v81/"
mutsig_run_script = "/home/andreas/bioinfo/core/mutsig/MutSigCV_1.4/run_MutSigCV.sh"

#############################################################################################################
# Genome MuSiC
#
# This tool suite doesn't need any hardcoded paths as it will install into the system path. 
# If installed correctly, typing "genome music" anywhere on the terminal should show the availabel commands.
#
# Check out these links:
# Paper: http://genome.cshlp.org/content/22/8/1589.long
# Homepage: http://gmt.genome.wustl.edu/genome-shipit/genome-music/0.4/index.html
# Installation https://www.biostars.org/p/62793/