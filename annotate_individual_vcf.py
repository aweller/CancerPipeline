import SampleAnnotation as SA
import sys
import logging

logging.basicConfig(level=logging.DEBUG,
                format='%(asctime)s %(name)-12s %(levelname)worker user-8s %(message)s',
                datefmt='%d-%m-%y %H:%M',
                #filename= "pipeline_log.txt",
                stream= sys.stdout)
                #filemode='w')

infile = sys.argv[1]

sample = SA.SampleAnnotation(infile, target_folder= "", run_tools=True)
sample.print_rows()