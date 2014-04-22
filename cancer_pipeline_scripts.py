import AnnotationRowParsers as AP
import SampleAnnotation as SA


def annotate_raw_vcfs(raw_vcf_folder):
    """
    Perform the initial Annovar/SNPeff annotation of unfiltered vcfs
    """
    
    raw_vcfs = [x for x in os.listdir(raw_vcf_folder)
                      if x.endswith(".vcf") and
                      "QC" not in x]

    for filename in raw_vcfs:
        #print "Starting", filename
        sample = SA.SampleAnnotation(filename, run_tools=True)
    SA.annotate_all_samples_as_one()