# import
try:
    import os, sys
    # if os.path.exists("/home/medhat/source/maf/alignio-maf"):
    #     sys.path.insert(1, "/home/medhat/source/maf/alignio-maf")
    # else:
    #     sys.path.insert(1, "/data/array1/medhat/artifical_read_project/genomes_assessment/alignio-maf")
    from RunAssessment import RunAssessment
    import argparse
except ImportError:
    print "oops, the import didn't work"

# code

parser = argparse.ArgumentParser(description='reading corrected file data and assessment the results')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.01')
parser.add_argument("-m", "--mafFile", help="A maf file that contains the alignment of original and artificial data",
                    type=str, required=True)
parser.add_argument("-c", "--correctedFile", help="A file that contains the corrected reads",
                    type=str, required=True)
parser.add_argument("-l", "--logFile", help="A PBcR log file that contains the corrected PacBio reads",
                    type=str, required=False)
parser.add_argument("-o", "--outputFile", help="A file that contains assessment of the corrected PacBio reads",
                    type=str, required=True)
parser.add_argument("-t", "--numberOfProcess", help="Number of running process default is 1",
                    type=int, required=False)
parser.add_argument("-d", "--dir", help="directory where work will be done",
                    type=str, required=True)
parser.add_argument("-s", "--software", help="software used for the corrected reads",
                    type=str, required=True)
parser.add_argument("-a", "--align", help="alignment output_type to be done",
                    type=str, required=True)


args = parser.parse_args()

maf = args.mafFile
corrected_file = args.correctedFile
log_map = args.logFile
if not log_map:
    log_map = None
result = args.outputFile
number_of_process = args.numberOfProcess
if not number_of_process:
    number_of_process = 1
working_dir = args.dir
software = args.software
align = args.align


if __name__ == "__main__":
    run_assessment = RunAssessment(corrected_fasta=corrected_file, maf_file=maf, corrected_file_type=software,
                                   assess_result_file=result, log_file=log_map, work_dir=working_dir,
                                   number_of_threads=number_of_process, align=align)
    run_assessment.correct_assess()
