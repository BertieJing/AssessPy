# import
from RunAssessment import RunAssessment
import argparse

parser = argparse.ArgumentParser(description='reading maf file data and assessment the alignment between reads')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.01')

parser.add_argument("-m", "--mafFile", help="A maf file that contains the alignment of original and artificial data",
                    type=str, required=True)
parser.add_argument("-t", "--numberOfProcess", help="Number of running process default is 1",
                    type=int, required=False)
parser.add_argument("-d", "--dir", help="directory where work will be done",
                    type=str, required=True)
parser.add_argument("-a", "--align", help="alignment output_type to be done",
                    type=str, required=True)
parser.add_argument("-o", "--outputFile", help="A file that contains assessment of the corrected PacBio reads",
                    type=str, required=True)

args = parser.parse_args()

maf = args.mafFile
result = args.outputFile
number_of_process = args.numberOfProcess
if not number_of_process:
    number_of_process = 1
working_dir = args.dir
align = args.align

if __name__ == '__main__':
    run_assess = RunAssessment.assess_maf(maf=maf, threads=number_of_process, directory=working_dir, output=result,
                                          align=align)

