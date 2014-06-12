import os
import subprocess
import re
from ToolConfig import vcf2maf_script, snpeff_jar, snpeff_config_file
import logging

def run_vcf2maf(input_vcf):
    
    if "/" in input_vcf: # samplename is a full path
        f = input_vcf.split("/")
        target_folder = "/".join(f[:-1]) + "/"
        input_vcf = f[-1]
    else:
        target_folder = "."
    
    output_maf = input_vcf[:-4] + ".maf"
    
    if os.path.exists(target_folder + output_maf):
        logging.debug( output_maf, "already present" )
    
    else:
        convert_cmd = "perl %s --input-vcf %s/%s --output-maf %s/%s " % (vcf2maf_script, target_folder, input_vcf, target_folder, output_maf)
        snpeff_cmd =  """--snpeff-cmd 'java -Xmx2g -jar %s \
                      hg19 -t -hgvs -c %s' """ % (snpeff_jar, snpeff_config_file) 
        
        logging.debug( convert_cmd+snpeff_cmd )
        subprocess.call(convert_cmd+snpeff_cmd, shell = True, stdout = open("log_out.txt", "wa"), stderr = open("log_err.txt", "wa"))
        
        try:
            os.remove(target_folder + input_vcf.strip(".vcf") +  ".anno.vcf")
        except:
            pass
        
        ###############################################
        # fix the TUMOR and NORMAL in the output
        # maf samples need to be of format sample.vcf (Tumor) and sample_N.vcf (Normal)
        
        sample = input_vcf.strip(".vcf").strip("_T").strip("_N")
        sample = re.sub("_v[1-9]", "", sample) # THIS NEEDS TO BE REMOVE AFTER THE FULL BAM NAMES ARE IMPORTED
        logging.warning(  "Removing _v[0-9] from sample names in module 'automate_vcf2maf.py'!" )

        tumor_name = sample
        normal_name = sample + "_N"
    
        sed_cmd = "sed -i 's/TUMOR/%s/g' %s/%s" % (tumor_name, target_folder, output_maf)
        subprocess.call(sed_cmd, shell = True)
    
        sed_cmd = "sed -i 's/NORMAL/%s/g' %s/%s" % (normal_name, target_folder, output_maf)
        subprocess.call(sed_cmd, shell = True)

def unite_mafs_in_current_folder():
    
    ###############################################
    # create united output file
    
    head_cmd = "cat *.maf | grep NCBI_Build | sort | uniq > tmp"
    cat_cmd = "cat *.maf | grep -v NCBI_Build >> tmp"
    mv_cmd = "mv tmp all_samples.maf"
    
    subprocess.call(head_cmd, shell=True)
    subprocess.call(cat_cmd, shell=True)
    subprocess.call(mv_cmd, shell=True)

def unite_mafs(samples, outname):
    
    if os.path.exists(outname):
        os.remove(outname)
    
    header = set([open(filename).readlines()[1] for filename in samples])
    assert len(header) == 1
    header = str(list(header)[0])

    all_rows =  [open(filename).readlines() for filename in samples]
    all_rows = list( set( sum(all_rows, []) )) # flatten list of lists    
    maf_rows = [x for x in all_rows if "NCBI_Build" not in x and not x.startswith("#")]
    
    with open(outname, "w") as out:
        out.write(str(header))
        for row in maf_rows:
            out.write(row)
        
def main():
    
    import sys
    import multiprocessing as mp

    cpus = int(sys.argv[1])
    
    to_run = [x for x in os.listdir(os.curdir) if x.endswith("vcf") and not "anno" in x and not "snpeff" in x]

    if cpus == 1:
        for vcf in to_run:
            logging.debug( "Converting %s to maf" % vcf)
            run_vcf2maf(vcf)
            
    else:
        pool = mp.Pool(cpus)
        for vcf in to_run:
            pool.apply_async(run_vcf2maf, args = [vcf,])
        pool.close()
        pool.join()        
        
    unite_mafs()

##############################################################################################

if __name__ == '__main__':
    #main()
    import sys
    
    infile = sys.argv[1]
    run_vcf2maf(infile)