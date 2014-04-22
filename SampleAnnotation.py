import os
import pprint
import sys
import subprocess
import AnnotationRowParsers as ARP
import automate_annovar as Autoanno
import automate_snpeff as Autosnp 
import collections

# TODO
# include Annovar parser that corrects for the imprecise chromposes

class SampleAnnotation():
    """
    Contains all rows of all annotation files for a given sample
    """
    
    def __init__(self, sample, target_folder = None, run_tools = False):
        
        if not target_folder:
            target_folder = ""
        
        if "/" in sample: # samplename is a full path
            f = sample.split("/")
            target_folder = "/".join(f[:-1]) + "/"
            sample = f[-1]
        
        self.sample = sample.split(".")[0]
        self.target_folder = target_folder
        
        self.vcf_rows = {}
        self.tsv_rows = {}
        self.anno_rows = {}
        
        ######################################################
        # initialize VCF
        vcf_filename = self.sample + ".vcf" 
        snpeff_filename = self.sample.split(".")[0] + ".snpeff.vcf"
        
        if not os.path.isfile(target_folder + snpeff_filename) and run_tools:
            #print "Cant find the Snpeff file, starting Snpeff."
            Autosnp.run_snpeff(target_folder, sample)
        
        if os.path.isfile(target_folder + snpeff_filename):
            vcf_filename = snpeff_filename
        
        ######################################################
        # initialize TSV
        tsv_filename_full = self.sample + "_full.tsv"
        tsv_filename_short = self.sample + ".tsv"
        
        if os.path.isfile(target_folder + tsv_filename_full):
            tsv_filename = tsv_filename_full
        elif os.path.isfile(target_folder + tsv_filename_short):
            tsv_filename = tsv_filename_short
        else:
            tsv_filename = "EMPTY"
        
        ######################################################
        # initialize ANNOVAR
        
        annovar_filename = self.sample + ".hg19_multianno.txt"
        if not os.path.isfile(target_folder + annovar_filename) and run_tools:
            #print "Cant find the Annovar file, starting Annovar."
            Autoanno.run_annovar(target_folder, sample)                    

        ######################################################
        # parse VCF
             
        with open(target_folder + vcf_filename) as vcf_file:
            for row in vcf_file:
                if "genotype" in row or row.startswith("#"): continue
                var = ARP.VCFrow(row)
                self.vcf_rows[var.chrompos] = var
#         print "parsed", target_folder + vcf_filename
        
        ######################################################
        # parse TSV
        
        if os.path.isfile(target_folder + tsv_filename):
            with open(target_folder + tsv_filename) as tsv_file:
                for row in tsv_file:
                    if "genotype" in row: continue
                    var = ARP.TSVrow(row)
                    self.tsv_rows[var.chrompos] = var
#             print "parsed", target_folder + tsv_filename
        
        ######################################################
        # parse ANNOVAR
        
        if os.path.isfile(target_folder + annovar_filename):        
            with open(target_folder + annovar_filename) as anno_file:
                for row in anno_file:
                    if "Start" in row: 
                        annovar_header = row.strip().split()
                    else:
                        var = ARP.ANNOVARrow(row, header = annovar_header)   
                        self.anno_rows[var.chrompos] = var
#             print "parsed", target_folder + annovar_filename

        ######################################################
        # remove intermediate files

        os.remove(target_folder + self.sample + ".hg19_multianno.txt")
        os.remove(target_folder + self.sample + ".snpeff.vcf")
        os.remove(target_folder + self.sample + ".annoin")
    
    def print_rows(self):
        vcf_keys, tsv_keys, anno_keys, all_keys = self._get_all_keys()
        all_chromposes = self.vcf_rows.keys() 
        all_chromposes.sort(cmp = sort_chromposes)
        
        output_filename = "%s_annotated.tsv" % self.sample
        output = open(self.target_folder + output_filename, "w")
        
        header_row = "\t".join(all_keys) + "\n"
        header_row = "sample\t" + header_row.replace("chrompos", "chrom\tpos")
        output.write(header_row)
        
        for chrompos in all_chromposes:
            vcf_row = self.vcf_rows[chrompos]
            vcf_values = [vcf_row.values.get(key, "-") for key in vcf_keys]
            
            tsv_values, anno_values = [], []
            
            tsv_row = self.tsv_rows.get(chrompos, None)
            if tsv_row:
                tsv_values = [tsv_row.values.get(key, "-") for key in tsv_keys]
            
            anno_row = self.anno_rows.get(chrompos)
            if anno_row:
                anno_values = [anno_row.values.get(key, "-") for key in anno_keys]
            
            result_values = [self.sample,] + vcf_values + tsv_values + anno_values
            result_row = "\t".join([str(x) for x in result_values]) + "\n"
            output.write(result_row)
        
        output.close()

    def _get_all_keys(self):
        vcf_keys = ["chrompos", "ref", "alt"]
        for var in self.vcf_rows.values():
            for key in var.values.keys():
                if key not in vcf_keys:
                    vcf_keys.append(key)
            
        tsv_keys = []
        for var in self.tsv_rows.values():
            for key in var.values.keys():
                if key not in tsv_keys:
                    tsv_keys.append(key)

        anno_keys = []
        for var in self.anno_rows.values():
            for key in var.values.keys():
                if key not in anno_keys:
                    anno_keys.append(key)
                    
        all_keys = vcf_keys + tsv_keys + anno_keys
        return vcf_keys, tsv_keys, anno_keys, all_keys
    
    #def _get_all_keys(self):
    #    core_keys = ["chrompos", "ref", "alt"]
    #    
    #    vcf_keys = [var.values.keys() for var in self.vcf_rows.values()] 
    #    vcf_keys = list( set( sum(vcf_keys, []) )) # flatten list of lists
    #    #vcf_keys.sort()
    #    vcf_keys = core_keys + [x for x in vcf_keys if x not in core_keys]
    #    
    #    tsv_keys = [var.values.keys() for var in self.tsv_rows.values()] 
    #    tsv_keys = list( set( sum(tsv_keys, []) )) # flatten list of lists
    #    #tsv_keys.sort()
    #
    #    anno_keys = [var.values.keys() for var in self.anno_rows.values()] 
    #    anno_keys = list( set( sum(anno_keys, []) )) # flatten list of lists
    #    #anno_keys.sort()
    #    
    #    all_keys = vcf_keys + tsv_keys + anno_keys
    #    return vcf_keys, tsv_keys, anno_keys, all_keys

def annotate_all_samples_as_one(filenames=None, folder=None, outfile=None):
        
    if not filenames and folder:
        filenames = [x for x in os.listdir(folder) if x.endswith(".vcf") and "snpeff" not in x and "anno" not in x]    
    elif not filenames and not folder:
        print "Please specify either a is of samples or a target folder."
    
    sample_no = len(filenames)
    
    target_folder = folder
    if not target_folder:
        target_folder = "./"
    
    if "/" in outfile: # samplename is a full path
        f = outfile.split("/")
        outfolder = "/".join(f[:-1]) + "/"
        outfile = f[-1]
    else:
        outfolder = "" 
    
    if outfile:
        final_outname= outfile
    else:
        final_outname= "all_samples_annotated.tsv" 
    
    # create a unified vcf of all variants in all samples
    
    all_rows =  [open(filename).readlines() for filename in filenames]
    all_rows = list( set( sum(all_rows, []) )) # flatten list of lists
    all_variant_rows = [row for row in all_rows if row[0] != "#"]
    #all_variant_rows.sort(cmp = sort_chromposes) #should not be necessary
    
    outname = outfolder + "all_samples.vcf" 
    with open(outname, "w") as out:
        for row in all_variant_rows:
            out.write(row)
            
    # annotate the unified vcf. this is kind of a waste as we already have these annotations, 
    # but cleaner than somehow passing the keys of all files around
    
    sample = SampleAnnotation(outname, run_tools = True)
    sample.print_rows()
    
    # recreate the united file, as the other file just uses "all_samples" as the samplename
    
    unite_cmd = "head -n 1 %s/all_samples_annotated.tsv > %s/%s && cat %s/*T_annotated.tsv | grep -v sample >> %s/%s" % (outfolder, outfolder, final_outname, outfolder, outfolder, final_outname)
    print unite_cmd
    subprocess.call(unite_cmd, shell=True)

    os.remove(outfolder + "all_samples.vcf")
    os.remove(outfolder + "all_samples_annotated.tsv")
    
def sort_chromposes(cp1, cp2):
    c1 = cp1.split("\t")[0].lstrip("chr")
    p1 = int(cp1.split("\t")[1])

    c2 = cp2.split("\t")[0].lstrip("chr")
    p2 = int(cp2.split("\t")[1])
    
    if c1 != c2:
        if c1 not in [str(x) for x in range(25)]:
            c1 = 98
        if c2 not in [str(x) for x in range(25)]:
            c2 = 99
            
        c1, c2 = int(c1), int(c2)
        if c1 > c2:
            result = 1
        elif c1 < c2:
            result = -1
    else:
        if p1 > p2:
            result = 1
        elif p1 < p2:
            return -1
        else:
            result = 1
    return result
    