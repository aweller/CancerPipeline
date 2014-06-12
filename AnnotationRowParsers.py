# This script contains the classes needed to parse VCFs, Annovar output and SNPeff output files

import sys
import re
import operator
import scipy.stats
import pprint
import collections
import logging

##############################################################################################################################################################################

class VCFrow():

    """
    Each instance corresponds to one variant site (= one row in a .vcf file)
    chr2    141771058    .    T    A    204.56    PASS    AO=95;DP=228;FAO=95;FRO=228;FR=.;FRO=133;FSAF=50;FSAR=45;FSRF=73;FSRR=60;FWDB=0.0166299;HRUN=4;LEN=1;MLLD=364.269;RBI=0.0288687;REFB=0;REVB=0.0235977;RO=133;SAF=50;SAR=45;SRF=71;SRR=62;SSEN=0;SSEP=0;STB=0.513213;SXB=0.0524288;TYPE=snp;VARB=0.0180844;OID=.;OPOS=141771058;OREF=T;OALT=A;OMAPALT=A    GT:GQ:DP:FRO:RO:FRO:AO:FAO:SAR:SAF:SRF:SRR:FSAR:FSAF:FSRF:FSRR    0/1:99:228:228:133:133:95:95:45:50:71:62:45:50:73:60
    'chr1\t115256528\t.\tT\t.\t100.0\tPASS\tHS;genes=NRAS;omim=164790;cosmic=585, 586, 587;dbsnp=rs121913255\tGT:GQ:GL:DP:FRO:AD:APSD:AST:ABQV\t0/0:90.22:-0.0,-999.0,-999.0:1437:1365:1363:726:34:27
    """

    def __init__(self, row):
        
        row = row.rstrip("\n")
        
        self.row = row
        self.fields = row.split("\t")
        
        self.values = collections.OrderedDict()
        self.chrompos = None
        
        if not row[0] == "#":
        
            self.chrom = self.fields[0]
            self.pos = int(self.fields[1])
            self.chrompos = self.chrom +"\t"+ str(self.pos)
    
            self.values["chrompos"] = self.chrompos
            self.values["ref"] = self.fields[3]
            self.values["alt"] = self.fields[4]
            
            self.values["called"] = self.fields[6]
        
            # info field
            self.info = {}
            self._fullinfo = self.fields[7]
            infofields = self._fullinfo.split(";")[1:] # because the first field indicates a hotspot
            for field in infofields:
                if "=" in field:
                    self.values[field.split("=")[0]]= field.split("=")[1] 
                    
            # format field
            self._formatkeys = self.fields[8]
            self._formatvalues = self.fields[9]
            self.format = {}
            for key,value in zip(self._formatkeys.split(":"), self._formatvalues.split(":")):
                
                if key == "FAO":
                    values = [int(x) for x in value.split(",")]
                    value = max(values) # as FAO can be a list for each allele
                self.values[key] = value
                
            # potential SNPeff field
            
            if self.values.get("EFF"):
                all_snpeffs =  self.values["EFF"].split(",")        
                all_snpeffs.sort(cmp=sort_snpeff_transcripts)
                snpeff = all_snpeffs[-1]
                
                effect = snpeff.split("(")[0]
                snpeff_fields = snpeff.split("(")[1].strip(")").split("|")
                
                self.values["Effect"] = effect
                
                snpeff_headers = ['Effect_Impact ', ' Functional_Class ', ' Codon_Change ', ' Amino_Acid_Change', ' Amino_Acid_Length ', ' Gene_Name ', ' Transcript_BioType ', ' Gene_Coding ', ' Transcript_ID ', ' Exon_Rank  ', ' Genotype_Number']
                snpeff_headers = [x.strip() for x in snpeff_headers]
                
                for key, value in zip(snpeff_headers, snpeff_fields):
                    self.values[key] = value
                
                self.values.pop("EFF")
                


def sort_snpeff_transcripts(t1, t2):
    eff_dict = {"HIGH":4, "MODERATE":3, "LOW":2, "MODIFIER":1}

    t1 = t1.split("(")[1].split("|")[0]
    t2 = t2.split("(")[1].split("|")[0] 
    
    eff1 = eff_dict[t1]
    eff2 = eff_dict[t2]

    if eff1 > eff2:
        return 1
    elif eff2 > eff1:
        return -1
    else:
        return 1

  
##############################################################################################################################################################################

class TSVrow():
    """
    Each instance corresponds to one row of an Ion Reporter tsv file
    1       231560226       A       A,C     EGLN1:  NM_022051.2:    utr_5:  :       :       :::     genes=EGLN1:dbsnp=rs1361383:omim=606425:pfam.acc=PF01
753:pfam.clan=TRASH:pfam.name=MYND finger:pfam.clan_id=CL0175:pfam.id=zf-MYND:pfam.panel=pfam:pfam.acc=PF03171:pfam.clan=Cupin:pfam.name=2OG-Fe(II) o
xygenase superfamily:pfam.clan_id=CL0029:pfam.id=2OG-FeII_Oxy:pfam.panel=pfam:phylop.sc=0.09
    """
    
    def __init__(self, row):

        row = row.rstrip("\n")
        
        self.row = row
        self.fields = row.split("\t")
        self.values = collections.OrderedDict()
        
        if not row[0] == "#":
            self.chrom = "chr" + self.fields[0]
            self.pos = int(self.fields[1])
            self.chrompos = self.chrom +"\t"+ str(self.pos)
    
            self.values["chrompos"] = self.chrompos
            self.values["ref"] = self.fields[2]
            self.values["alt"] = self.fields[3]
            
            #chr    pos     ref     genotype        genes   transcripts     location        function        proteinChange   funcScores      custom
            
            self.values["genes"] = self.fields[4]
            self.values["transcripts"] = self.fields[5]
            self.values["location"] = self.fields[6]
            self.values["function"] = self.fields[7]
            self.values["proteinChange"] = self.fields[8]
            self.values["funcScores"] = self.fields[9]
            
            self.custom = self.fields[10]
            
            self.customfields = [x.split(";") for x in self.custom.split(":")]
            self.customfields = list( set( sum(self.customfields, []) )) # flatten list of lists
            for item in self.customfields:
                if not "=" in item: continue
                key, value = item.split("=")
                self.values[key] = value

##############################################################################################################################################################################

class ANNOVARrow():
    
    def __init__(self, row, header = None, filename = None):

        row = row.rstrip("\n")
        
        self.row = row
        self.fields = row.split("\t")
        self.values = collections.OrderedDict()
              
        self.chrom = self.fields[0]
        self.pos = int(self.fields[1])
        self.chrompos = self.chrom +"\t"+ str(self.pos)        
        
        if not header and not filename:
            logging.critical( "Can't parse Annovar output without header row (or filename to parse it from)." )
            logging.critical(  "Exiting..." )
            sys.exit()
        
        if not header: # have to fetch header
            with open(filename) as handle:
                for row in handle:
                    if row[0] != "#":
                        header = row.strip().split("\t")
                        break                
        
        if not row[0] == "#":
            for key, value in zip(header, self.fields):
                self.values[key] = value
        

        
        