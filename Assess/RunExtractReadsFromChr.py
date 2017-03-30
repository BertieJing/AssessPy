""" Script to extract uncorrected  and corrected reads coordinate from chromosome
"""
import pickle
import argparse
from ExtractReadsFromChr import ExtractReadsFromChr

parser = argparse.ArgumentParser(
    description='Script to extract uncorrected  and corrected reads coordinate from chromosome')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.01')
parser.add_argument("-m", "--mafFile", help="A maf file that contains the alignment of original and artificial data",
                    type=str, required=True)
parser.add_argument("-c", "--correctedReads", help="trimmed splitted proovread file",
                    type=str, required=True)
parser.add_argument("-d", "--dir", help="directory where work will be done",
                    type=str, required=True)

args = parser.parse_args()
maf = args.mafFile
corrected = args.correctedReads
working_dir = args.dir

if __name__ == "__main__":
    a = ExtractReadsFromChr(maf_file=maf, corrected_file=corrected, working_dir=working_dir)
    corrected, uncorrected = a.get_corrected_uncorrected_coordinate()
    pickle.dump(corrected, open(working_dir+"/corrected_dir_coordinate.p", "rb"))
    pickle.dump(uncorrected, open(working_dir+"/uncorrected_dir_coordinate.p", "rb"))
    print('Done...')