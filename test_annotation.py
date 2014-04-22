import AnnotationRowParsers as AP
import SampleAnnotation as SA
import sys

infile = sys.argv[1]

sample = SA.SampleAnnotation(infile, target_folder= "", run_tools=True)
sample.print_rows()